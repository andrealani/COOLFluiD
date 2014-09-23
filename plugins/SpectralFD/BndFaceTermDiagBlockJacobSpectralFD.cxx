#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/BndFaceTermDiagBlockJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BndFaceTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule >
  BndFaceTermDiagBlockJacobSpectralFDProvider("BndFaceTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

BndFaceTermDiagBlockJacobSpectralFD::BndFaceTermDiagBlockJacobSpectralFD(const std::string& name) :
  BndFaceTermRHSJacobSpectralFD(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_currDiagMatrix(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

BndFaceTermDiagBlockJacobSpectralFD::~BndFaceTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermDiagBlockJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  BndFaceTermRHSJacobSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermDiagBlockJacobSpectralFD::executeOnTrs()
{
  CFAUTOTRACE;

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

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

          // ADD TOTAL UPDATE TO UPDATE COEFFICIENTS
          m_resUpdates += m_diffResUpdates;
          m_updateCoeffContr += m_diffUpdateCoeffContr;
          addUpdateCoeffContribution();

          // COMPUTE JACOBIAN CONTRIBUTION OF FACE TERMS
          m_currDiagMatrix = &diagBlockJacobMatr[m_intCell->getID()];
          computeJacobBndFaceTerm();
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermDiagBlockJacobSpectralFD::computeJacobBndFaceTerm()
{
  // get number of solution points in the cell
  const CFuint nbrSolPnts = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // loop over the states in the internal cell to perturb the states
  CFuint pertResUpdIdx = 0;
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert, ++pertResUpdIdx)
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
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable and its gradient in the flux points
      restoreFacePhysVarsAndGrad(iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes <= m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  BndFaceTermRHSSpectralFD::setup();

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // maximum number of solution points in a cell
  const CFuint maxNbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  const CFuint nbrRes = maxNbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrRes);
  m_derivResUpdates.resize(nbrRes);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  BndFaceTermRHSSpectralFD::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    BndFaceTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = BndFaceTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD
