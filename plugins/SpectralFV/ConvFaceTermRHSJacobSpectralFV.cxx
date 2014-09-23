#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFV/ConvFaceTermRHSJacobSpectralFV.hh"
#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SpectralFVElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ConvFaceTermRHSJacobSpectralFV, SpectralFVMethodData, SpectralFVModule> ConvFaceTermRHSJacobSpectralFVProvider("ConvFaceTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermRHSJacobSpectralFV::ConvFaceTermRHSJacobSpectralFV(const std::string& name) :
  ConvFaceTermRHSSpectralFV(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_derivResUpdates()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermRHSJacobSpectralFV::~ConvFaceTermRHSJacobSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFV::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFV::configure ( Config::ConfigArgs& args )
{
  ConvFaceTermRHSSpectralFV::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFV::execute()
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
        m_faceTermComputer->reconstructFluxPntsStates(*m_states[LEFT],*m_states[RIGHT]);
        m_faceTermComputer->reconstructFaceAvgState  (*m_states[LEFT],*m_states[RIGHT]);
      }

      // if one of the neighbouring cells is parallel updatable,
      // compute the face term
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // compute the face terms and the wave speed updates
        m_faceTermComputer
                      ->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,
                                                              m_waveSpeedUpd[LEFT ],
                                                              m_waveSpeedUpd[RIGHT]);

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

void ConvFaceTermRHSJacobSpectralFV::computeBothJacobsConvFaceTerm()
{
  // number of subfaces
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  CFuint cvIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrCVs = m_states[iSide]->size();
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV, ++cvIdx)
    {
      acc.setRowColIndex(cvIdx,(*m_states[iSide])[iCV]->getLocalID());
    }
  }

  // number of CVs in left cell
  const CFuint nbrLeftCVs = m_states[LEFT]->size();

  // loop over left and right cell
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // loop over the states/CVs in the left and right cell to perturb the states
    const CFuint nbrCVs = m_states[iSide]->size();
    for (CFuint iCVPert = 0; iCVPert < nbrCVs; ++iCVPert)
    {
      // dereference state
      State& pertState = *(*m_states[iSide])[iCVPert];

      // loop over the variables in the state
      for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
      {
        // perturb physical variable in state
        m_numJacob->perturb(iEqPert,pertState[iEqPert]);

        // backup and reconstruct physical variable in the flux points
        m_faceTermComputer->backupAndReconstructPhysVar(iSide,iEqPert,*m_states[iSide]);

        // compute the perturbed face term
        m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);

        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;

        // add the derivative of the residual updates to the accumulator
        const CFuint pertSideTerm = iSide*nbrLeftCVs;
        // right cell
        CFuint resUpdIdx = 0;
        for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx +=m_nbrEqs)
        {
          CFuint cvIdx = (*m_cvCVConnSVFace)[m_orient][iSubFace][RIGHT];
          acc.addValues(cvIdx+nbrLeftCVs,iCVPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
        }
        // left cell
        m_derivResUpdates *= -1.0; // different sign for left cell
        resUpdIdx = 0;
        for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx +=m_nbrEqs)
        {
          CFuint cvIdx = (*m_cvCVConnSVFace)[m_orient][iSubFace][LEFT ];
          acc.addValues(cvIdx,iCVPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
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

void ConvFaceTermRHSJacobSpectralFV::computeOneJacobConvFaceTerm(const CFuint side)
{
  // factor depending on whether dealing with left or right cell
  const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(side+1));

  // number of subfaces
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  CFuint cvIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrCVs = m_states[iSide]->size();
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV, ++cvIdx)
    {
      acc.setRowColIndex(cvIdx,(*m_states[iSide])[iCV]->getLocalID());
    }
  }

  // number of CVs in left cell
  const CFuint nbrLeftCVs = m_states[LEFT]->size();

  // term depending on the side
  const CFuint sideTerm = side*nbrLeftCVs;

  // loop over left and right cell
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // loop over the states/CVs to perturb the states
    const CFuint nbrCVs = m_states[iSide]->size();
    for (CFuint iCVPert = 0; iCVPert < nbrCVs; ++iCVPert)
    {
      // dereference states
      State& pertState = *(*m_states[iSide])[iCVPert];

      // loop over the variables in the state
      for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
      {
        // perturb physical variable in state
        m_numJacob->perturb(iEqPert,pertState[iEqPert]);

        // backup and reconstruct physical variable in the flux points
        m_faceTermComputer->backupAndReconstructPhysVar(iSide,iEqPert,*m_states[iSide]);

        // compute the perturbed face term
        m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);

        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;

        // add the derivative of the residual updates to the accumulator
        m_derivResUpdates *= sideFactor;
        const CFuint pertSideTerm = iSide*nbrLeftCVs;
        CFuint resUpdIdx = 0;
        for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx +=m_nbrEqs)
        {
          CFuint cvIdx = (*m_cvCVConnSVFace)[m_orient][iSubFace][side];
          acc.addValues(cvIdx+sideTerm,iCVPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
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

void ConvFaceTermRHSJacobSpectralFV::setup()
{
  CFAUTOTRACE;

  ConvFaceTermRHSSpectralFV::setup();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  cf_assert(svLocalData.size() > 0);

  // get the number of CVs in a cell
  const CFuint nbrCVs = svLocalData[0]->getNbrOfCVs();

  // get SV face - CV connectivity
  SafePtr< vector< vector< CFuint > > > svFaceCVConn = svLocalData[0]->getExtSVFaceCVConn();

  // maximum number of CVs at SV face
  CFuint maxNbrCVsAtSVFace = 0;
  for (CFuint iFace = 0; iFace < svFaceCVConn->size(); ++iFace)
  {
    const CFuint nbrCVsAtSVFace = (*svFaceCVConn)[iFace].size();
    maxNbrCVsAtSVFace = maxNbrCVsAtSVFace > nbrCVsAtSVFace ? maxNbrCVsAtSVFace : nbrCVsAtSVFace;
  }

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(2*nbrCVs,2*nbrCVs,m_nbrEqs));

  // resize variables
  const CFuint nbrFaceFluxes = maxNbrCVsAtSVFace*m_nbrEqs;
  m_pertResUpdates .resize(nbrFaceFluxes);
  m_derivResUpdates.resize(nbrFaceFluxes);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSJacobSpectralFV::unsetup()
{
  CFAUTOTRACE;

  ConvFaceTermRHSSpectralFV::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
