#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/DiffFaceTermDiagBlockJacobSpectralFD.hh"
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

MethodCommandProvider<DiffFaceTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> DiffFaceTermDiagBlockJacobSpectralFDProvider("DiffFaceTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

DiffFaceTermDiagBlockJacobSpectralFD::DiffFaceTermDiagBlockJacobSpectralFD(const std::string& name) :
  DiffFaceTermRHSJacobSpectralFD(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_currDiagMatrix()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffFaceTermDiagBlockJacobSpectralFD::~DiffFaceTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermDiagBlockJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

  // set the data needed to compute the face terms;
  setFaceTermData();

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

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

  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCBL = m_cellBuilders[LEFT ]->getDataGE();
  geoDataCBL.trs = cells;
  CellToFaceGEBuilder::GeoData& geoDataCBR = m_cellBuilders[RIGHT]->getDataGE();
  geoDataCBR.trs = cells;

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

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_face->getNeighborGeo(LEFT )->getStates();
      m_states[RIGHT] = m_face->getNeighborGeo(RIGHT)->getStates();
      const bool lParUpdatable = (*m_states[LEFT ])[0]->isParUpdatable();
      const bool rParUpdatable = (*m_states[RIGHT])[0]->isParUpdatable();

      // if one of the neighbouring cells is parallel updatable,
      // compute the face terms
      if (lParUpdatable || rParUpdatable)
      {
        // build the neighbouring cells
        const CFuint cellIDL = m_face->getNeighborGeo(LEFT )->getID();
        geoDataCBL.idx = cellIDL;
        m_cells[LEFT ] = m_cellBuilders[LEFT ]->buildGE();
        const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
        geoDataCBR.idx = cellIDR;
        m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();

        // set the gradients
        setGradients();

        // set the current face and compute the face data in the face term computer
        m_faceTermComputer->setCurrentFace(m_face);
        m_faceTermComputer->computeFaceData();
        m_faceTermComputer->computeNeighbourCellData();

        // reconstruct the states
        m_faceTermComputer->reconstructFluxPntsStates(m_states);

        // reconstruct the gradients
        m_faceTermComputer->reconstructFluxPntsGradients(m_grads);

        // compute the face terms and the update coefficient contributions
        m_faceTermComputer
            ->computeDiffFaceTermAndUpdateCoefContributions(m_resUpdates,m_updateCoefContr);

        // add update coefficient contributions
        addUpdateCoeffContributions();

        // get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current boundary faces
        setOtherFacesLocalIdxs();

        // get the neigbouring states of the other faces
        setFaceNeighbourStates();

        // get the neigbouring gradients of the other faces
        setFaceNeighbourGradients();

        // set the neighbouring cells data
        setCellsData();

        // compute auxiliary term for the perturbed gradient reconstruction from current face
        computeCellGradsMinusFaceTerm();

        // compute the diffusive face term contribution to the jacobian
        if (lParUpdatable)
        {
          m_currDiagMatrix = &diagBlockJacobMatr[m_cells[LEFT ]->getID()];
          computeOneJacobDiffFaceTerm(LEFT );
        }
        if (rParUpdatable)
        {
          m_currDiagMatrix = &diagBlockJacobMatr[m_cells[RIGHT]->getID()];
          computeOneJacobDiffFaceTerm(RIGHT);
        }

        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
      }

      // release the face
      m_faceBuilder->releaseGE();
    }
  }
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermDiagBlockJacobSpectralFD::computeOneJacobDiffFaceTerm(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // variable for the other side
  const CFuint otherSide = side == LEFT ? RIGHT : LEFT;

  // loop over the states to perturb the states
  const CFuint nbrSolPnts = m_states[side]->size();
  CFuint pertResUpdIdx = 0;
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_states[side])[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert, ++pertResUpdIdx)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);

      // backup and reconstruct physical variable in the face flux points
      m_faceTermComputer->backupAndReconstructPhysVar(side,iEqPert,*m_states[side]);

      // backup and reconstruct physical variable in the other face flux points
      // and in the cell flux points
      backupAndReconstructOtherFacesAndCellPhysVars(side,iEqPert);

      // compute the perturbed gradients in the current cell
      computePerturbedGradients(side);

      // compute the perturbed gradients in the other cell
      computePertGradsFromFaceTerm(otherSide);

      // reconstruct the face gradients
      m_faceTermComputer->reconstructFluxPntsGradients(m_pertGrads);

      // compute the perturbed face term
      m_faceTermComputer->computeDiffFaceTerm(m_pertResUpdates);

      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates[side],m_resUpdates[side],m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to diagonal block matrix
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points (face term computer)
      m_faceTermComputer->restorePhysVar(side,iEqPert);

      // restore physical variable in the other face flux points
      // and in the cell flux points
      restoreOtherFacesAndCellPhysVars(side,iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes <= m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class of parent class, in order not to call the lss
  DiffFaceTermRHSSpectralFD::setup();

  // get CellToFaceGeBuilders
  m_cellBuilders.resize(2);
  m_cellBuilders[LEFT ] = getMethodData().getCellBuilder();
  m_cellBuilders[RIGHT] = getMethodData().getSecondCellBuilder();

  // get some additional data for cell building
  m_isFaceOnBoundary.resize(2);
  m_nghbrCellSide   .resize(2);
  m_currCellSide    .resize(2);
  m_faceOrients     .resize(2);
  m_faceBCIdx       .resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_isFaceOnBoundary[iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getIsFaceOnBoundary();
    m_nghbrCellSide   [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getNeighbrCellSide ();
    m_currCellSide    [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getCurrentCellSide ();
    m_faceOrients     [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getFaceOrient      ();
    m_faceBCIdx       [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getFaceBCIdx       ();
  }

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get number of equations and dimensionality
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // get the number of solution points in a cell
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // get the number of faces in a cell
  const CFuint nbrFaces = sdLocalData[0]->getNbrCellFaces();
  const CFuint nbrFacesM1 = nbrFaces - 1;

  // get face term computers and boundary face term computers
  m_faceTermComputers   .resize(2);
  m_bndFaceTermComputers.resize(2);
  CFuint computerIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_faceTermComputers   [iSide].resize(nbrFacesM1);
    m_bndFaceTermComputers[iSide].resize(nbrFacesM1);
    for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace, ++computerIdx)
    {
      m_faceTermComputers   [iSide][iFace] = getMethodData().getAdditionalFaceTermComputer   (computerIdx);
      m_bndFaceTermComputers[iSide][iFace] = getMethodData().getAdditionalBndFaceTermComputer(computerIdx);
    }
  }

  // get the volume term computers
  m_volTermComputers.resize(2);
  m_volTermComputers[LEFT ] = getMethodData().getVolTermComputer();
  m_volTermComputers[RIGHT] = getMethodData().getSecondVolTermComputer();

  // resize m_faces
  m_faces.resize(2);

  // resize m_otherFaceLocalIdxs
  m_otherFaceLocalIdxs.resize(2);
  m_otherFaceLocalIdxs[LEFT ].resize(nbrFacesM1);
  m_otherFaceLocalIdxs[RIGHT].resize(nbrFacesM1);

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_faceNghbrStates[iSide].resize(nbrFacesM1);
    for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
    {
      m_faceNghbrStates[iSide][iFace].resize(2);
    }
  }

  // resize m_faceNghbrStates
  m_faceNghbrGrads.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_faceNghbrGrads[iSide].resize(nbrFacesM1);
    for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
    {
      m_faceNghbrGrads[iSide][iFace].resize(2);
    }
  }

  // resize variables
  const CFuint nbrCellResiduals = nbrSolPnts*m_nbrEqs;
  m_pertResUpdates.resize(2);
  m_pertResUpdates[LEFT ].resize(nbrCellResiduals);
  m_pertResUpdates[RIGHT].resize(nbrCellResiduals);
  m_derivResUpdates.resize(nbrCellResiduals);
  m_unpertCellDiffRes       .resize(2);
  m_unpertCellDiffRes[LEFT ].resize(nbrCellResiduals);
  m_unpertCellDiffRes[RIGHT].resize(nbrCellResiduals);
  m_pertCellDiffRes         .resize(nbrCellResiduals);
  m_derivCellDiffRes        .resize(nbrCellResiduals);

  // allocate memory for perturbed gradiens
  m_pertGrads.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_pertGrads[iSide].resize(nbrSolPnts);
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_pertGrads[iSide][iSol] = new vector<RealVector>(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        (*m_pertGrads[iSide][iSol])[iEq].resize(dim);
      }
    }
  }

  // resize m_cellGradsMinusFaceTerm
  m_cellGradsMinusFaceTerm.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_cellGradsMinusFaceTerm[iSide].resize(nbrSolPnts);
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_cellGradsMinusFaceTerm[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_cellGradsMinusFaceTerm[iSide][iSol][iEq].resize(dim);
      }
    }
  }

  // resize m_cellGradsMinusOtherFaceTerm
  m_cellGradsMinusOtherFaceTerm.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_cellGradsMinusOtherFaceTerm[iSide].resize(nbrFacesM1);
    for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
    {
      m_cellGradsMinusOtherFaceTerm[iSide][iFace].resize(nbrSolPnts);
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        m_cellGradsMinusOtherFaceTerm[iSide][iFace][iSol].resize(m_nbrEqs);
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_cellGradsMinusOtherFaceTerm[iSide][iFace][iSol][iEq].resize(dim);
        }
      }
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

  // resize neighbouring cells solution points Jacobian determinants
  m_solJacobDet.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_solJacobDet[iSide].resize(nbrSolPnts);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSide = 0; iSide < m_pertGrads.size(); ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_pertGrads[iSide].size(); ++iSol)
    {
      deletePtr(m_pertGrads[iSide][iSol]);
    }
    m_pertGrads[iSide].resize(0);
  }
  m_pertGrads.resize(0);

  // call unsetup of parent class of parent class, in order not to call the lss
  DiffFaceTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffFaceTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffFaceTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
