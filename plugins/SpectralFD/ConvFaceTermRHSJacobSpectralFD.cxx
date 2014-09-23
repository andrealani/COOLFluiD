#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "SpectralFD/ConvFaceTermRHSJacobSpectralFD.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ConvFaceTermRHSJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> ConvFaceTermRHSJacobSpectralFDProvider("ConvFaceTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermRHSJacobSpectralFD::ConvFaceTermRHSJacobSpectralFD(const std::string& name) :
  ConvFaceTermRHSSpectralFD(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_derivResUpdates()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermRHSJacobSpectralFD::~ConvFaceTermRHSJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  ConvFaceTermRHSSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

  // set the data needed to compute the face terms;
  setFaceTermData();

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;
  geoData.isBoundary = false;

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // set the orientation of the faces
    m_faceTermComputer->setFaceOrientation(m_orient);

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoData.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      // if one of the neighbouring cells is parallel updatable or the gradients have to be computed,
      // set face data and reconstruct states
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable() || hasDiffTerm)
      {
        // set the current face and compute the face data in the face term computer
        m_faceTermComputer->setCurrentFace(m_face);
        m_faceTermComputer->computeFaceData();

        // reconstruct the states
        m_faceTermComputer->reconstructFluxPntsStates(m_states);
      }

      // if one of the neighbouring cells is parallel updatable,
      // compute the face term
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // compute the face terms and the wave speed updates
        m_faceTermComputer
            ->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_waveSpeedUpd);

        // update the rhs
        updateRHS();

        // update the wave speeds
        updateWaveSpeed();
      }

      // compute the convective face term contribution to the jacobian
      if ((*m_states[LEFT ])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
      {
        computeBothJacobsConvFaceTerm();
      }
      else if ((*m_states[LEFT ])[0]->isParUpdatable())
      {
        computeOneJacobConvFaceTerm(LEFT );
      }
      else if ((*m_states[RIGHT])[0]->isParUpdatable())
      {
        computeOneJacobConvFaceTerm(RIGHT);
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradientFaceTerm();
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFD::computeBothJacobsConvFaceTerm()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // number of solution points in left cell
  const CFuint nbrLeftSolPnts = m_states[LEFT]->size();

  // set block row and column indices
  CFuint solIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[iSide])[iSol]->getLocalID());
    }
  }

  // loop over left and right cell
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // term depending on iSide
    const CFuint pertSideTerm = iSide*nbrLeftSolPnts;

    // loop over the states in the left and right cell to perturb the states
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
    {
      // dereference state
      State& pertState = *(*m_states[iSide])[iSolPert];

      // loop over the variables in the state
      for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
      {
        // perturb physical variable in state
        m_numJacob->perturb(iEqPert,pertState[iEqPert]);

        // backup and reconstruct physical variable in the flux points
        m_faceTermComputer->backupAndReconstructPhysVar(iSide,iEqPert,*m_states[iSide]);

        // compute the perturbed face term
        m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);

        // add contributions to the Jacobian
        for (CFuint iSide2 = 0; iSide2 < 2; ++iSide2)
        {
          // factor depending on iSide2
          const CFuint sideTerm = iSide2*nbrLeftSolPnts;

          // compute the finite difference derivative of the face term
          m_numJacob->computeDerivative(m_pertResUpdates[iSide2],m_resUpdates[iSide2],m_derivResUpdates);

          // multiply residual update derivatives with residual factor
          m_derivResUpdates *= resFactor;

          // add the derivative of the residual updates to the accumulator
          CFuint resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+sideTerm,iSolPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[iEqPert]);

        // restore physical variable in the flux points
        m_faceTermComputer->restorePhysVar(iSide,iEqPert);
      }
    }
  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFD::computeOneJacobConvFaceTerm(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  CFuint solIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[iSide])[iSol]->getLocalID());
    }
  }

  // number of solution points in left cell
  const CFuint nbrLeftSolPnts = m_states[LEFT]->size();

  // term depending on the side
  const CFuint sideTerm = side*nbrLeftSolPnts;

  // loop over left and right cell
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // term depending on iSide
    const CFuint pertSideTerm = iSide*nbrLeftSolPnts;

    // loop over the states to perturb the states
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
    {
      // dereference states
      State& pertState = *(*m_states[iSide])[iSolPert];

      // loop over the variables in the state
      for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
      {
        // perturb physical variable in state
        m_numJacob->perturb(iEqPert,pertState[iEqPert]);

        // backup and reconstruct physical variable in the flux points
        m_faceTermComputer->backupAndReconstructPhysVar(iSide,iEqPert,*m_states[iSide]);

        // compute the perturbed face term
        m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);

        // add contributions to the Jacobian
        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates[side],m_resUpdates[side],m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;

        // add the derivative of the residual updates to the accumulator
        CFuint resUpdIdx = 0;
        for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
        {
          acc.addValues(iSol+sideTerm,iSolPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[iEqPert]);

        // restore physical variable in the flux points
        m_faceTermComputer->restorePhysVar(iSide,iEqPert);
      }
    }
  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  ConvFaceTermRHSSpectralFD::setup();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // maximum number of solution points in a cell
  const CFuint maxNbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(2*maxNbrSolPnts,2*maxNbrSolPnts,m_nbrEqs));

  // resize variables
  const CFuint resSize = maxNbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(2);
  m_pertResUpdates [LEFT ].resize(resSize);
  m_pertResUpdates [RIGHT].resize(resSize);
  m_derivResUpdates.resize(resSize);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  ConvFaceTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
