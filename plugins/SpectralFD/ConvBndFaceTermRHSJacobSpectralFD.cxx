#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/ConvBndFaceTermRHSJacobSpectralFD.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndFaceTermRHSJacobSpectralFD, SpectralFDMethodData, SpectralFDModule >
  ConvBndFaceTermRHSJacobSpectralFDProvider("ConvBndFaceTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSJacobSpectralFD::ConvBndFaceTermRHSJacobSpectralFD(const std::string& name) :
  ConvBndFaceTermRHSSpectralFD(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_derivResUpdates()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSJacobSpectralFD::~ConvBndFaceTermRHSJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  ConvBndFaceTermRHSSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSJacobSpectralFD::executeOnTrs()
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

  // get current QuadFreeBCSpectralFD TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get bndFacesStartIdxs from SpectralFDMethodData
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

void ConvBndFaceTermRHSJacobSpectralFD::computeJacobConvBndFaceTerm()
{
  // get number of solution points in the cell
  const CFuint nbrSolPnts = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
  }

  // loop over the states in the internal cell to perturb the states
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iSolPert];

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

      // add contributions to the Jacobian
      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
      {
        acc.addValues(iSol,iSolPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
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

void ConvBndFaceTermRHSJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  ConvBndFaceTermRHSSpectralFD::setup();

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
  m_acc.reset(m_lss->createBlockAccumulator(maxNbrSolPnts,maxNbrSolPnts,m_nbrEqs));

  // resize variables
  const CFuint nbrRes = maxNbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrRes);
  m_derivResUpdates.resize(nbrRes);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  ConvBndFaceTermRHSSpectralFD::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD
