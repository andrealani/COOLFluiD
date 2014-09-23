#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/DiffBndFaceTermDiagBlockJacobSpectralFD.hh"
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

MethodCommandProvider< DiffBndFaceTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule >
  DiffBndFaceTermDiagBlockJacobSpectralFDProvider("DiffBndFaceTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

DiffBndFaceTermDiagBlockJacobSpectralFD::DiffBndFaceTermDiagBlockJacobSpectralFD(const std::string& name) :
  DiffBndFaceTermRHSJacobSpectralFD(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_currDiagMatrix(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndFaceTermDiagBlockJacobSpectralFD::~DiffBndFaceTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermDiagBlockJacobSpectralFD::executeOnTrs()
{
  CFAUTOTRACE;

  // set BCStateComputer in the boundary face term computer
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // set the data needed to compute the face terms;
  setFaceTermData();

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

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

  // get the geodata of the cell builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCB = m_cellBuilder->getDataGE();
  geoDataCB.trs = cellTrs;

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

        // get the states in the neighbouring cell
        m_cellStates = m_face->getNeighborGeo(0)->getStates();

        // if cell is parallel updatable, compute the boundary face term
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // build the neighbouring cell
          const CFuint cellID = m_face->getNeighborGeo(0)->getID();
          geoDataCB.idx = cellID;
          m_intCell = m_cellBuilder->buildGE();

          // set current diagonal block Jacobian matrix
          m_currDiagMatrix = &diagBlockJacobMatr[m_intCell->getID()];

          // set the gradients
          setGradients();

          // set the current face and compute the face data in the boundary face term computer
          m_bndFaceTermComputer->setCurrentFace(m_face);
          m_bndFaceTermComputer->computeFaceData();
          m_bndFaceTermComputer->computeNeighbourCellData();

          // reconstruct the states
          m_bndFaceTermComputer->reconstructFluxPntsStates(*m_cellStates);

          // reconstruct the gradients
          m_bndFaceTermComputer->reconstructFluxPntsGradients(m_cellGrads);

          // compute the face terms and the update coefficient contributions
          m_bndFaceTermComputer->computeDiffFaceTermAndUpdateCoefContributions(m_resUpdates,m_updateCoeffContr);

          // update the contribution to the update coefficient
          addUpdateCoeffContribution();

          // get all the faces neighbouring the cell
          m_faces = m_intCell->getNeighborGeos();

          // set the local indexes of the other faces than the current boundary faces
          setOtherFacesLocalIdxs();

          // get the neigbouring states of the other faces
          setFaceNeighbourStates();

          // set data for the neighbouring cell
          setCellData();

          // compute the contribution to the jacobian
          computeJacobDiffBndFaceTerm();

          // release the cell
          m_cellBuilder->releaseGE();
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermDiagBlockJacobSpectralFD::computeJacobDiffBndFaceTerm()
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

      // backup and reconstruct physical variable in the boundary face flux points
      // and reconstruct the ghost states
      m_bndFaceTermComputer->backupAndReconstructPhysVar(iEqPert,*m_cellStates);

      // backup and reconstruct physical variable in the other face flux points
      // and in the cell flux points
      backupAndReconstructOtherFacesAndCellPhysVars(iEqPert);

      // recompute the cell gradients
      computePerturbedGradients();

      // reconstruct the boundary face gradients
      m_bndFaceTermComputer->reconstructFluxPntsGradients(m_pertGrads);

      // compute the perturbed boundary face term
      m_bndFaceTermComputer->computeDiffFaceTerm(m_pertResUpdates);

      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to diagonal block matrix
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the boundary face flux points
      m_bndFaceTermComputer->restorePhysVar(iEqPert);

      // restore physical variable in the other face flux points
      // and in the cell flux points
      restoreOtherFacesAndCellPhysVars(iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes <= m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class of parent class, in order not to call the lss
  DiffBndFaceTermRHSSpectralFD::setup();

  // get CellToFaceGeBuilder
  m_cellBuilder      = getMethodData().getCellBuilder();
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSide    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSide     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrients      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdx        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  // get the volume term computer
  m_volTermComputer = getMethodData().getVolTermComputer();

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get number of equations and dimensionality
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points in a cell
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // get the number of faces in a cell
  const CFuint nbrFaces = sdLocalData[0]->getNbrCellFaces();
  const CFuint nbrFacesM1 = nbrFaces - 1;

  // get face term computers and boundary face term computers
  m_faceTermComputers   .resize(nbrFacesM1);
  m_bndFaceTermComputers.resize(nbrFacesM1);
  for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
  {
    m_faceTermComputers   [iFace] = getMethodData().getAdditionalFaceTermComputer   (iFace);
    m_bndFaceTermComputers[iFace] = getMethodData().getAdditionalBndFaceTermComputer(iFace);
  }

  // resize m_otherFaceLocalIdxs
  m_otherFaceLocalIdxs.resize(nbrFacesM1);

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(nbrFacesM1);
  for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

  // resize variables
  const CFuint nbrCellResiduals = nbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrCellResiduals);
  m_derivResUpdates.resize(nbrCellResiduals);

  // allocate memory for perturbed gradients
  m_pertGrads.resize(nbrSolPnts);
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_pertGrads[iSol] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[iSol])[iEq].resize(dim);
    }
  }

  // resize gradient updates
  m_gradUpdates.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_gradUpdates[iSide].resize(nbrSolPnts);
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_gradUpdates[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_gradUpdates[iSide][iSol][iEq].resize(dim);
      }
    }
  }

  // resize neighbouring cell solution point Jacobian determinants
  m_solJacobDet.resize(nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSol = 0; iSol < m_pertGrads.size(); ++iSol)
  {
    deletePtr(m_pertGrads[iSol]);
  }
  m_pertGrads.resize(0);

  // call unsetup of parent class of parent class, in order not to call the lss
  DiffBndFaceTermRHSSpectralFD::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffBndFaceTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffBndFaceTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD
