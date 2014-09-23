#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/ConvBndFaceTermRHSJacobSpectralFV.hh"
#include "SpectralFV/SpectralFV.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndFaceTermRHSJacobSpectralFV, SpectralFVMethodData, SpectralFVModule >
  ConvBndFaceTermRHSJacobSpectralFVProvider("ConvBndFaceTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSJacobSpectralFV::ConvBndFaceTermRHSJacobSpectralFV(const std::string& name) :
  ConvBndFaceTermRHSSpectralFV(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_derivResUpdates()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSJacobSpectralFV::~ConvBndFaceTermRHSJacobSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSJacobSpectralFV::configure ( Config::ConfigArgs& args )
{
  ConvBndFaceTermRHSSpectralFV::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSJacobSpectralFV::executeOnTrs()
{
  CFAUTOTRACE;

  // set BCStateComputer in the boundary face term computer
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // set the data needed to compute the face terms;
  setFaceTermData();

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCSpectralFV TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get bndFacesStartIdxs from SpectralFVMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // set the orientation of the faces
      m_bndFaceTermComputer->setFaceOrientation(m_orient);

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);

        // get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();

        // if cell is parallel updatable or the gradients have to be computed,
        // set face data and reconstruct states
        if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
        {
          // set the current face and compute the face data in the boundary face term computer
          m_bndFaceTermComputer->setCurrentFace(m_face);
          m_bndFaceTermComputer->computeFaceData();

          // reconstruct the states
          m_bndFaceTermComputer->reconstructFluxPntsStates(*m_cellStates);
          m_bndFaceTermComputer->reconstructFaceAvgState  (*m_cellStates);
        }

        // if cell is parallel updatable, compute the boundary face term
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // compute the face terms and the wave speed updates
          m_bndFaceTermComputer->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_waveSpeedUpd);

          // update the rhs
          updateRHS();

          // update the wave speeds in the neighbouring cell
          updateWaveSpeed();

          // compute the convective boundary face term contribution to the jacobian
          computeJacobConvBndFaceTerm();
        }

        // if there is a diffusive term, compute the gradients
        if (hasDiffTerm)
        {
          computeGradientFaceTerm();
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSJacobSpectralFV::computeJacobConvBndFaceTerm()
{
  // get number of CVs in the cells
  const CFuint nbrCVs = m_cellStates->size();

  // number of subfaces
  const CFuint nbrSubFaces = (*m_svFaceCVConn)[m_orient].size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    acc.setRowColIndex(iCV,(*m_cellStates)[iCV]->getLocalID());
  }

  // loop over the states/CVs in the internal cell to perturb the states
  for (CFuint iCVPert = 0; iCVPert < nbrCVs; ++iCVPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iCVPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);

      // backup and reconstruct physical variable in the flux points
      // and reconstruct the ghost states)
      m_bndFaceTermComputer->backupAndReconstructPhysVar(iEqPert,*m_cellStates);

      // compute the perturbed boundary face term
      m_bndFaceTermComputer->computeConvFaceTerm(m_pertResUpdates);

      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_resUpdates,m_pertResUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx += m_nbrEqs)
      {
        CFuint cvIdx = (*m_svFaceCVConn)[m_orient][iSubFace];
        acc.addValues(cvIdx,iCVPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points
      m_bndFaceTermComputer->restorePhysVar(iEqPert);
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

void ConvBndFaceTermRHSJacobSpectralFV::setup()
{
  CFAUTOTRACE;

  ConvBndFaceTermRHSSpectralFV::setup();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  cf_assert(svLocalData.size() > 0);

  // get SV face - CV connectivity
  SafePtr< vector< vector< CFuint > > > svFaceCVConn = svLocalData[0]->getExtSVFaceCVConn();

  // maximum number of CVs at SV face
  CFuint maxNbrCVsAtSVFace = 0;
  for (CFuint iFace = 0; iFace < svFaceCVConn->size(); ++iFace)
  {
    const CFuint nbrCVsAtSVFace = (*svFaceCVConn)[iFace].size();
    maxNbrCVsAtSVFace = maxNbrCVsAtSVFace > nbrCVsAtSVFace ? maxNbrCVsAtSVFace : nbrCVsAtSVFace;
  }

  // get the number of CVs in a cell
  const CFuint nbrCVs = svLocalData[0]->getNbrOfCVs();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(nbrCVs,nbrCVs,m_nbrEqs));

  // resize variables
  const CFuint nbrFaceFluxes = maxNbrCVsAtSVFace*m_nbrEqs;
  m_pertResUpdates .resize(nbrFaceFluxes);
  m_derivResUpdates.resize(nbrFaceFluxes);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSJacobSpectralFV::unsetup()
{
  CFAUTOTRACE;

  ConvBndFaceTermRHSSpectralFV::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFV

} // namespace COOLFluiD
