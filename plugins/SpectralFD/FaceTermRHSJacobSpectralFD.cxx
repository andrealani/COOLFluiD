#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "SpectralFD/FaceTermRHSJacobSpectralFD.hh"
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

MethodCommandProvider<FaceTermRHSJacobSpectralFD, SpectralFDMethodData, SpectralFDModule>
  FaceTermRHSJacobSpectralFDProvider("FaceTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

FaceTermRHSJacobSpectralFD::FaceTermRHSJacobSpectralFD(const std::string& name) :
  FaceTermRHSSpectralFD(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_derivResUpdates()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

FaceTermRHSJacobSpectralFD::~FaceTermRHSJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  FaceTermRHSSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

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

  // SET THE DATA NEEDED TO COMPUTE THE FACE TERMS
  setFaceTermComputerData();

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // SET THE ORIENTATION OF THE FACES
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

      // if one of the neighbouring cells is parallel updatable, compute the face terms
      const bool lParUpdatable = (*m_states[LEFT ])[0]->isParUpdatable();
      const bool rParUpdatable = (*m_states[RIGHT])[0]->isParUpdatable();
      if (lParUpdatable || rParUpdatable)
      {
        // COMPUTE DATA IN FACE TERM COMPUTER
        computeFaceTermData();

        // COMPUTE CONVECTIVE FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
        m_faceTermComputer->
            computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_updateCoefContr);

        // COMPUTE DIFFUSIVE FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
        m_faceTermComputer->
            computeDiffFaceTermAndUpdateCoefContributions(m_diffResUpdates,m_diffUpdateCoefContr);

        // ADD TOTAL UPDATE TO RESIDUAL AND UPDATE COEFFICIENTS
        m_resUpdates[LEFT ] += m_diffResUpdates[LEFT ];
        m_resUpdates[RIGHT] += m_diffResUpdates[RIGHT];
        addUpdatesToResidual();
        m_updateCoefContr[LEFT ] += m_diffUpdateCoefContr[LEFT ];
        m_updateCoefContr[RIGHT] += m_diffUpdateCoefContr[RIGHT];
        addUpdateCoeffContributions();

        // COMPUTE JACOBIAN CONTRIBUTION OF FACE TERMS
        if (lParUpdatable && rParUpdatable)
        {
          computeBothJacobsFaceTerm();
        }
        else if (lParUpdatable)
        {
          computeOneJacobFaceTerm(LEFT);
        }
        else if (rParUpdatable)
        {
          computeOneJacobFaceTerm(RIGHT);
        }
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
// CF_DEBUG_EXIT;
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSJacobSpectralFD::backupAndReconstructFacePhysVarsAndGrad(const CFuint side, const CFuint iVar)
{
  // physical variable
  m_faceTermComputer->backupAndReconstructPhysVar(side,iVar,*m_states[side]);

  // physical variable gradient
  m_faceTermComputer->backupAndReconstrFluxPntsSolPolyGrad(side,iVar,*m_states[side]);
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSJacobSpectralFD::restoreFacePhysVarsAndGrad(const CFuint side, const CFuint iVar)
{
  // physical variable
  m_faceTermComputer->restorePhysVar(side,iVar);

  // physical variable gradient
  m_faceTermComputer->restorePhysVarGrad(side,iVar);
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSJacobSpectralFD::computeBothJacobsFaceTerm()
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

        // backup and reconstruct physical variable and its gradient in the flux points
        backupAndReconstructFacePhysVarsAndGrad(iSide,iEqPert);

        // DEBUG
//         computeFaceTermData();

        // compute the perturbed face term
        m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);
        m_faceTermComputer->computeDiffFaceTerm(m_diffResUpdates);

        // add contributions to the Jacobian
        for (CFuint iSide2 = 0; iSide2 < 2; ++iSide2)
        {
          // factor depending on iSide2
          const CFuint sideTerm = iSide2*nbrLeftSolPnts;

          // compute the finite difference derivative of the face term
          m_pertResUpdates[iSide2] += m_diffResUpdates[iSide2];
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

        // restore physical variable and its gradient in the flux points
        restoreFacePhysVarsAndGrad(iSide,iEqPert);
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

void FaceTermRHSJacobSpectralFD::computeOneJacobFaceTerm(const CFuint side)
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

        // backup and reconstruct physical variable and its gradient in the flux points
        backupAndReconstructFacePhysVarsAndGrad(iSide,iEqPert);

        // compute the perturbed face term
        m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);
        m_faceTermComputer->computeDiffFaceTerm(m_diffResUpdates);

        // add contributions to the Jacobian
        // compute the finite difference derivative of the face term
        m_pertResUpdates[side] += m_diffResUpdates[side];
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

        // restore physical variable and its gradient in the flux points
        restoreFacePhysVarsAndGrad(iSide,iEqPert);
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

void FaceTermRHSJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  FaceTermRHSSpectralFD::setup();

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

void FaceTermRHSJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of parent class
  FaceTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
