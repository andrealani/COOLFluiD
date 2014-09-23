#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/BndFaceTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BndFaceTermRHSJacobSpectralFD, SpectralFDMethodData, SpectralFDModule >
  BndFaceTermRHSJacobSpectralFDProvider("BndFaceTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

BndFaceTermRHSJacobSpectralFD::BndFaceTermRHSJacobSpectralFD(const std::string& name) :
  BndFaceTermRHSSpectralFD(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_derivResUpdates()
{
}

//////////////////////////////////////////////////////////////////////////////

BndFaceTermRHSJacobSpectralFD::~BndFaceTermRHSJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  BndFaceTermRHSSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSJacobSpectralFD::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current TRS
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

  // SET THE DATA NEEDED TO COMPUTE THE FACE TERMS
  setFaceTermComputerData();

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // SET THE ORIENTATION OF THE FACES
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

        // if cell is parallel updatable, compute the boundary face term
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceTermData();

          // COMPUTE CONVECTIVE BOUNDARY FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
          m_bndFaceTermComputer->
              computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_updateCoeffContr);

          // COMPUTE DIFFUSIVE BOUNDARY FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
          m_bndFaceTermComputer->
              computeDiffFaceTermAndUpdateCoefContributions(m_diffResUpdates,m_diffUpdateCoeffContr);

          // ADD TOTAL UPDATE TO RESIDUAL AND UPDATE COEFFICIENTS
          m_resUpdates += m_diffResUpdates;
          addUpdatesToResidual();
          m_updateCoeffContr += m_diffUpdateCoeffContr;
          addUpdateCoeffContribution();

          // COMPUTE JACOBIAN CONTRIBUTION OF FACE TERMS
          computeJacobBndFaceTerm();
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSJacobSpectralFD::backupAndReconstructFacePhysVarsAndGrad(const CFuint iVar)
{
  // physical variable
  m_bndFaceTermComputer->backupAndReconstructPhysVar(iVar,*m_cellStates);

  // physical variable gradient
  m_bndFaceTermComputer->backupAndReconstrFluxPntsSolPolyGrad(iVar,*m_cellStates);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSJacobSpectralFD::restoreFacePhysVarsAndGrad(const CFuint iVar)
{
  // physical variable
  m_bndFaceTermComputer->restorePhysVar(iVar);

  // physical variable gradient
  m_bndFaceTermComputer->restorePhysVarGrad(iVar);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSJacobSpectralFD::computeJacobBndFaceTerm()
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

      // backup and reconstruct physical variable and its gradient in the flux points
      // and reconstruct the ghost states and ghost gradients
      backupAndReconstructFacePhysVarsAndGrad(iEqPert);

      // compute the perturbed boundary face term
      m_bndFaceTermComputer->computeConvFaceTerm(m_pertResUpdates);
      m_bndFaceTermComputer->computeDiffFaceTerm(m_diffResUpdates);

      // add contributions to the Jacobian
      // compute the finite difference derivative of the face term
      m_pertResUpdates += m_diffResUpdates;
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

      // restore physical variable and its gradient in the flux points
      restoreFacePhysVarsAndGrad(iEqPert);
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

void BndFaceTermRHSJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  BndFaceTermRHSSpectralFD::setup();

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

void BndFaceTermRHSJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  BndFaceTermRHSSpectralFD::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD
