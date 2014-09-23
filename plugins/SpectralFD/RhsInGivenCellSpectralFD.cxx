#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/RhsInGivenCellSpectralFD.hh"
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

MethodCommandProvider<RhsInGivenCellSpectralFD, SpectralFDMethodData, SpectralFDModule>
    RhsInGivenCellSpectralFDProvider("RhsInGivenCell");

//////////////////////////////////////////////////////////////////////////////

RhsInGivenCellSpectralFD::RhsInGivenCellSpectralFD(const string& name) :
  SpectralFDMethodCom(name),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_updateCoeff("updateCoeff"),
  socket_statesSetIdx("statesSetIdx"),
  m_cellBuilder(CFNULL),
  m_nghbrCellBuilder(CFNULL),
  m_volTermComputer(CFNULL),
  m_volTermComputerNghbrCell(CFNULL),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_faceTermComputersNghbrCell(),
  m_bndFaceTermComputersNghbrCell(),
  m_bcStateComputers(),
  m_cell(CFNULL),
  m_nghbrCell(CFNULL),
  m_cellStates(CFNULL),
  m_nghbrCellStates(CFNULL),
  m_faceNghbrStates(),
  m_isFaceOnBoundary(CFNULL),
  m_nghbrCellSide(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_faceNghbrStatesNghbrCell(CFNULL),
  m_isFaceOnBoundaryNghbrCell(CFNULL),
  m_nghbrCellSideNghbrCell(CFNULL),
  m_currCellSideNghbrCell(CFNULL),
  m_faceOrientsNghbrCell(CFNULL),
  m_faceBCIdxNghbrCell(CFNULL),
  m_resUpdates(),
  m_faceResUpdates(),
  m_updateCoefUpd(),
  m_updateCoefFaceUpd(),
  m_nbrEqs(),
  m_dim(),
  m_iElemType(),
  m_currIter(),
  m_computeUpdateCoef(),
  m_gradUpdates(),
  m_currCellGrads(),
  m_nghbCellGrads(),
  m_solJacobDet(),
  m_solJacobDetNghbrCell(),
  m_solPntsLocalCoords(CFNULL),
  m_otherFaceLocalIdxs()
{
}

//////////////////////////////////////////////////////////////////////////////

RhsInGivenCellSpectralFD::~RhsInGivenCellSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::execute()
{
  // update current iteration number
  const CFuint currIter = SubSystemStatusStack::getActive()->getNbIter();
  if (m_currIter != currIter)
  // if we are at a new time step, then the diagonal block Jacobian matrices could
  // have been computed, and the data in the volume and face term computers should be reset.
  {
    // update m_currIter
    m_currIter = currIter;

    // set the volume and face term computers' data
    setVolumeAndFaceTermComputersData();

    // resize the residual updates if necessary
    resizeResUpdates();
  }

  // get index of cell for which to compute the rhs
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  const CFuint cellIdx = statesSetIdx[0];
  m_computeUpdateCoef = statesSetIdx[1];

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the cell builder and set the TRS and the cell index
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  geoData.idx = cellIdx;

  // build the cell
  m_cell = m_cellBuilder->buildGE();

  // get the states in this cell
  m_cellStates = m_cell->getStates();

  // if cell is parallel updatable, compute the rhs
  if ((*m_cellStates)[0]->isParUpdatable())
  {
    /// @note KVDA: for now, it is assumed that there is only one element type in the mesh.
    /// Could pass the element type from LUSGSMethod, in the same socket as the statesSetIdx.
    const CFuint currElemType = 0;
    if (m_iElemType != currElemType)
    {
      m_iElemType = currElemType;

      // set the volume and face term computers' data
      setVolumeAndFaceTermComputersData();

      // resize the residual updates if necessary
      resizeResUpdates();
    }

    // CLEAR THE CURRENT CELL RESIDUAL
    clearResidual();

    // COMPUTE THE VOLUME TERM
    setVolumeTermDataAndComputeConvVolumeTerm();

    // COMPUTE THE FACE TERMS
    setFaceTermDataAndComputeConvFaceTermsAndUpdateCoefs();

    // ADD DIFFUSIVE RESIDUAL IF NECESSARY
    if (getMethodData().hasDiffTerm())
    {
      // compute the gradients in the current cell and the diffusive volume terms
      computeCurrentCellGradientsAndDiffVolumeTerm();

      // compute the diffusive face terms
      computeDiffFaceTermsAndUpdateCoefs();
    }

    // MULTIPLY RHS WITH RESIDUAL FACTOR
    multplyRHSWithResFactor();
  }

  //release the GeometricEntity
  m_cellBuilder->releaseGE();
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::setVolumeAndFaceTermComputersData()
{
  // set the volume term data in the volume term computer
  m_volTermComputer->setVolumeTermData(m_iElemType);

  // set the data in the face term and boundary face term computers
  for (CFuint iFace = 0; iFace < m_faceTermComputers.size(); ++iFace)
  {
    m_faceTermComputers[iFace]->setFaceTermData();
  }
  for (CFuint iFace = 0; iFace < m_bndFaceTermComputers.size(); ++iFace)
  {
    m_bndFaceTermComputers[iFace]->setFaceTermData();
  }

  // set data for diffusive residuals
  if (getMethodData().hasDiffTerm())
  {
    // get the local spectral FD data
    vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

    // get solution point local coordinates
    m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::resizeResUpdates()
{
  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // resize m_resUpdates
  const CFuint nbrRes = nbrSolPnts*m_nbrEqs;
  m_resUpdates.resize(nbrRes);
  m_faceResUpdates[LEFT ].resize(nbrRes);
  m_faceResUpdates[RIGHT].resize(nbrRes);

  // resize m_gradUpdates
  if (getMethodData().hasDiffTerm())
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_gradUpdates[iSide].resize(nbrSolPnts);
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        m_gradUpdates[iSide][iSol].resize(m_nbrEqs);
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_gradUpdates[iSide][iSol][iEq].resize(m_dim);
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::clearResidual()
{
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();
  rhsCurrStatesSet = 0.0;

  if (m_computeUpdateCoef)
  {
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    const CFuint nbrSolPnts = m_cellStates->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
      updateCoeff[stateID] = 0.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::addUpdatesToResidual()
{
  // get the datahandle of the current cell rhs
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

  // update current cell rhs
  const CFuint nbrRes = m_resUpdates.size();
  cf_assert(nbrRes <= rhsCurrStatesSet.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    rhsCurrStatesSet[iRes] += m_resUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::addUpdateToUpdateCoef()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // add wave speed to current cell update coefficients
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
    updateCoeff[stateID] += m_updateCoefUpd;
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::setVolumeTermDataAndComputeConvVolumeTerm()
{
  // set the current cell and compute the cell data in the volume term computer
  m_volTermComputer->setCurrentCell(m_cell);
  m_volTermComputer->computeCellData();

  // reconstruct the solution in the flux points
  m_volTermComputer->reconstructStates(*m_cellStates);

  // compute the convective volume term
  m_volTermComputer->computeCellConvVolumeTerm(m_resUpdates);

  // add volume term to the current cell rhs
  addUpdatesToResidual();
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::setFaceTermDataAndComputeConvFaceTermsAndUpdateCoefs()
{
  const vector< Framework::GeometricEntity* >* faces = m_cell->getNeighborGeos();

  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    // set face neighbour states
    const CFuint nbrFaceNeighbours = (*faces)[iFace]->nbNeighborGeos();
    for (CFuint iSide = 0; iSide < nbrFaceNeighbours; ++iSide)
    {
      GeometricEntity* cell = (*faces)[iFace]->getNeighborGeo(iSide);
      m_faceNghbrStates[iFace][iSide] = cell->getStates();
    }

    if ((*m_isFaceOnBoundary)[iFace])
    {
      // set BC state computer
      m_bndFaceTermComputers[iFace]->setBcStateComputer((*m_bcStateComputers)[(*m_faceBCIdx)[iFace]]);

      // set the orientation of the face (should be set here, might be changed in the command that computes the face terms)
      m_bndFaceTermComputers[iFace]->setFaceOrientation(iFace);

      // set the face in the boundary face term computer
      m_bndFaceTermComputers[iFace]->setCurrentFace((*faces)[iFace]);

      // compute the face data in the boundary face term computer
      m_bndFaceTermComputers[iFace]->computeFaceData();

      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->reconstructFluxPntsStates(*m_cellStates);

      // compute the convective face term
      if (m_computeUpdateCoef)
      {
        m_bndFaceTermComputers[iFace]
          ->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_updateCoefUpd);
      }
      else
      {
        m_bndFaceTermComputers[iFace]->computeConvFaceTerm(m_resUpdates);
      }
    }
    else
    {
      // set the orientation of the face
      m_faceTermComputers[iFace]->setFaceOrientation((*m_faceOrients)[iFace]);

      // set the face in the face term computer
      m_faceTermComputers[iFace]->setCurrentFace((*faces)[iFace]);

      // compute the face data in the face term computer
      m_faceTermComputers[iFace]->computeFaceData();

      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->reconstructFluxPntsStates(m_faceNghbrStates[iFace]);

      // get current cell side
      const CFuint currCellSide = (*m_currCellSide)[iFace];

      // compute the convective face term
      if (m_computeUpdateCoef)
      {
        m_faceTermComputers[iFace]
            ->computeConvFaceTermAndWaveSpeedUpdates(m_faceResUpdates,m_updateCoefFaceUpd);
        m_updateCoefUpd = m_updateCoefFaceUpd[currCellSide];
      }
      else
      {
        m_faceTermComputers[iFace]->computeConvFaceTerm(m_faceResUpdates);
      }

      // copy the updates for the current cell
      m_resUpdates = m_faceResUpdates[currCellSide];
    }

    // add face term to the current cell rhs
    addUpdatesToResidual();

    // add update to current cell update coefficient
    if (m_computeUpdateCoef)
    {
      addUpdateToUpdateCoef();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::multplyRHSWithResFactor()
{
  // get factor for the residual
  const CFreal resFactor = getMethodData().getResFactor();

  // get the datahandle of the current cell rhs
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

  // update current cell rhs
  const CFuint nbrRes = m_resUpdates.size();
  cf_assert(nbrRes <= rhsCurrStatesSet.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    rhsCurrStatesSet[iRes] *= resFactor;
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::computeCurrentCellGradientsAndDiffVolumeTerm()
{
  // compute volume contribution to current cell gradients
  m_volTermComputer->computeGradientVolumeTerm(m_gradUpdates[0]);
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_currCellGrads[iSol])[iEq] = m_gradUpdates[0][iSol][iEq];
    }
  }

  // compute face contributions to the gradients
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // compute the boundary face contribution to the gradients
      m_bndFaceTermComputers[iFace]->computeGradientFaceTerm(m_gradUpdates[0]);

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_currCellGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
        }
      }
    }
    else
    {
      // compute the face contribution to the gradients
      m_faceTermComputers[iFace]->computeGradientFaceTerm(m_gradUpdates);

      // add the contribution to the gradients
      const CFuint side = (*m_currCellSide)[iFace];
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_currCellGrads[iSol])[iEq] += m_gradUpdates[side][iSol][iEq];
        }
      }
    }
  }

  // compute solution points Jacobian determinants
  m_solJacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // divide by solution point Jacobian determinant
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDet[iSol];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_currCellGrads[iSol])[iEq] *= invJacobDet;
    }
  }

  // reconstruct the gradients in the flux points
  m_volTermComputer->reconstructGradients(m_currCellGrads);

  // compute the diffusive volume term
  m_volTermComputer->computeCellDiffVolumeTerm(m_resUpdates);

  // add volume term to the current cell rhs
  addUpdatesToResidual();
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::computeDiffFaceTermsAndUpdateCoefs()
{
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // reconstruct the gradients
      m_bndFaceTermComputers[iFace]->reconstructFluxPntsGradients(m_currCellGrads);

      // compute the face terms and the update coefficient contributions
      if (m_computeUpdateCoef)
      {
        // compute the neighbouring cell data in the boundary face term computer (volume, needed for diffusive update coefficient)
        m_bndFaceTermComputers[iFace]->computeNeighbourCellData();

        m_bndFaceTermComputers[iFace]
            ->computeDiffFaceTermAndUpdateCoefContributions(m_resUpdates,m_updateCoefUpd);
      }
      else
      {
        m_bndFaceTermComputers[iFace]->computeDiffFaceTerm(m_resUpdates);
      }
    }
    else
    {
      // build neighbouring cell
      buildNeighbourCell(iFace);

      // set local indexes of the other faces
      setOtherFacesLocalIdxs(iFace);

      // set the neighbouring states to faces of the current neighbour cell
      setFaceNeighbourStates();

      // set the data in the neigbouring cell volume term computer and face term computers
      setNeighbourCellData();

      // compute the gradients in the neighbouring cell
      computeNeighbourCellGradients(iFace);

      // get current cell and other cell sides
      const CFuint currCellSide  = (*m_currCellSide)[iFace];
      const CFuint nghbCellSide = currCellSide == LEFT ? RIGHT : LEFT;

      // reconstruct the gradients
      m_faceTermComputers[iFace]->reconstructFluxPntsGradients(currCellSide,m_currCellGrads);
      m_faceTermComputers[iFace]->reconstructFluxPntsGradients(nghbCellSide,m_nghbCellGrads);

      // compute the face terms and the update coefficient contributions
      if (m_computeUpdateCoef)
      {
        // compute the neighbouring cell data in the face term computer (volume, needed for diffusive update coefficient)
        m_faceTermComputers[iFace]->computeNeighbourCellData();

        m_faceTermComputers[iFace]->computeDiffFaceTermAndUpdateCoefContributions(m_faceResUpdates,m_updateCoefFaceUpd);
        m_updateCoefUpd = m_updateCoefFaceUpd[currCellSide];
      }
      else
      {
        m_faceTermComputers[iFace]->computeDiffFaceTerm(m_faceResUpdates);
      }

      // copy the updates for the current cell
      m_resUpdates = m_faceResUpdates[currCellSide];

      //release the current neighbouring cell
      m_nghbrCellBuilder->releaseGE();
    }

    // add face term to the current cell rhs
    addUpdatesToResidual();

    // add update to current cell update coefficient
    if (m_computeUpdateCoef)
    {
      addUpdateToUpdateCoef();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::buildNeighbourCell(const CFuint currFaceIdx)
{
  // get neighbour cell ID
  const vector< Framework::GeometricEntity* >* faces = m_cell->getNeighborGeos();
  cf_assert((*faces)[currFaceIdx]->nbNeighborGeos() == 2);
  cf_assert((*faces)[currFaceIdx]->getNeighborGeo(LEFT )->getID() == m_cell->getID() ||
                   (*faces)[currFaceIdx]->getNeighborGeo(RIGHT)->getID() == m_cell->getID());
  const CFuint nghbrCellID =
          (*faces)[currFaceIdx]->getNeighborGeo(LEFT )->getID() == m_cell->getID() ?
          (*faces)[currFaceIdx]->getNeighborGeo(RIGHT)->getID() :
          (*faces)[currFaceIdx]->getNeighborGeo(LEFT )->getID() ;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the neighbour cell builder and set the TRS and the cell index
  CellToFaceGEBuilder::GeoData& geoData = m_nghbrCellBuilder->getDataGE();
  geoData.trs = cells;
  geoData.idx = nghbrCellID;

  // build the cell
  m_nghbrCell = m_nghbrCellBuilder->buildGE();

  // get the states in this cell
  m_nghbrCellStates = m_nghbrCell->getStates();
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::setOtherFacesLocalIdxs(const CFuint currFaceIdx)
{
  // get face ID of current face
  const vector< Framework::GeometricEntity* >* faces = m_cell->getNeighborGeos();
  const CFuint currFaceID = (*faces)[currFaceIdx]->getID();

  // find neighbour cell neighbouring faces != current face
  const CFuint nbrFaces = m_nghbrCell->nbNeighborGeos();
  const vector< Framework::GeometricEntity* >* nghbrCellFaces = m_nghbrCell->getNeighborGeos();
  cf_assert(nghbrCellFaces->size() == nbrFaces);
  CFuint iFace = 0;
  for (CFuint faceIdx = 0; faceIdx < nbrFaces; ++faceIdx)
  {
    if ((*nghbrCellFaces)[faceIdx]->getID() != currFaceID)
    {
      m_otherFaceLocalIdxs[iFace] = faceIdx;
      ++iFace;
    }
  }
  cf_assert(iFace <= m_otherFaceLocalIdxs.size());
  cf_assert(iFace == nbrFaces-1);
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::setFaceNeighbourStates()
{
  // get neighbour states of other faces
  const vector< Framework::GeometricEntity* >* nghbrCellFaces = m_nghbrCell->getNeighborGeos();
  const CFuint nbrOtherFaces = nghbrCellFaces->size() - 1;
  cf_assert(nbrOtherFaces <= m_otherFaceLocalIdxs.size());
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
      // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];

      // get neigbouring states
    const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundaryNghbrCell)[faceIdx] ? 1 : 2;
    for (CFuint iSide = 0; iSide < nbrFaceNeighbours; ++iSide)
    {
      GeometricEntity* cell = (*nghbrCellFaces)[faceIdx]->getNeighborGeo(iSide);
      m_faceNghbrStatesNghbrCell[iFace][iSide] = cell->getStates();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::setNeighbourCellData()
{
  // compute solution points Jacobian determinants
  /// @note KVDA: here it is assumed that the neighbouring cells all have the same type as the current cell
  /// by using the same values for m_solPntsLocalCoords, and m_iElemType for setting the volume term data
  m_solJacobDetNghbrCell =
      m_nghbrCell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // set the neighbouring cell and compute neighbouring cell data in the volume term computer
  m_volTermComputerNghbrCell->setVolumeTermData(m_iElemType);
  m_volTermComputerNghbrCell->setCurrentCell(m_nghbrCell);
  m_volTermComputerNghbrCell->computeCellData();
  m_volTermComputerNghbrCell->reconstructStates(*m_nghbrCellStates);

  // set orientation or boundary condition state computers of the other faces
  const vector< Framework::GeometricEntity* >* nghbrCellFaces = m_nghbrCell->getNeighborGeos();
  const CFuint nbrOtherFaces = nghbrCellFaces->size() - 1;
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];

    // set orientation or boundary condition
    if ((*m_isFaceOnBoundaryNghbrCell)[faceIdx])
    {
      // set data for computation of the face terms
      m_bndFaceTermComputersNghbrCell[iFace]->setFaceTermData();

      // set BC state computer
      m_bndFaceTermComputersNghbrCell[iFace]
          ->setBcStateComputer((*m_bcStateComputers)[(*m_faceBCIdxNghbrCell)[faceIdx]]);

      // set the orientation of the face
      m_bndFaceTermComputersNghbrCell[iFace]->setFaceOrientation(faceIdx);

      // set the face in the boundary face term computer
      m_bndFaceTermComputersNghbrCell[iFace]->setCurrentFace((*nghbrCellFaces)[faceIdx]);

      // compute the face data in the boundary face term computer
      m_bndFaceTermComputersNghbrCell[iFace]->computeFaceData();

      // compute the neighbouring cell data in the boundary face term computer
//       m_bndFaceTermComputersNghbrCell[iFace]->computeNeighbourCellData();

      // reconstruct the states in the face flux points
      m_bndFaceTermComputersNghbrCell[iFace]->reconstructFluxPntsStates(*m_nghbrCellStates);
    }
    else
    {
      // set data for computation of the face terms
      m_faceTermComputersNghbrCell[iFace]->setFaceTermData();

      // set the orientation of the face
      m_faceTermComputersNghbrCell[iFace]->setFaceOrientation((*m_faceOrientsNghbrCell)[faceIdx]);

      // set the face in the face term computer
      m_faceTermComputersNghbrCell[iFace]->setCurrentFace((*nghbrCellFaces)[faceIdx]);

      // compute the face data in the face term computer
      m_faceTermComputersNghbrCell[iFace]->computeFaceData();

      // compute the neighbouring cell data in the face term computer
//       m_faceTermComputersNghbrCell[iFace]->computeNeighbourCellData();

      // reconstruct the states in the face flux points
      m_faceTermComputersNghbrCell[iFace]->
          reconstructFluxPntsStates(m_faceNghbrStatesNghbrCell[iFace]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::computeNeighbourCellGradients(const CFuint currFaceIdx)
{
  // COMPUTE VOLUME CONTRIBUTION TO NEIGHBOUR CELL GRADIENTS
  m_volTermComputerNghbrCell->computeGradientVolumeTerm(m_gradUpdates[0]);
  const CFuint nbrSolPnts = m_nghbrCellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_nghbCellGrads[iSol])[iEq] = m_gradUpdates[0][iSol][iEq];
    }
  }

  // COMPUTE FACE CONTRIBUTIONS (OF THE OTHER FACES) TO THE GRADIENTS
  const vector< Framework::GeometricEntity* >* nghbrCellFaces = m_nghbrCell->getNeighborGeos();
  const CFuint nbrOtherFaces = nghbrCellFaces->size() - 1;
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];
    if ((*m_isFaceOnBoundaryNghbrCell)[faceIdx])
    {
      // compute the boundary face contribution to the gradients
      m_bndFaceTermComputersNghbrCell[iFace]->computeGradientFaceTerm(m_gradUpdates[0]);

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_nghbCellGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
        }
      }
    }
    else
    {
      // compute the face contribution to the gradients
      m_faceTermComputersNghbrCell[iFace]->computeGradientFaceTerm(m_gradUpdates);

      // add the contribution to the gradients
      const CFuint side = (*m_currCellSideNghbrCell)[faceIdx];
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_nghbCellGrads[iSol])[iEq] += m_gradUpdates[side][iSol][iEq];
        }
      }
    }
  }

  // COMPUTE FACE CONTRIBUTION (OF THE CURRENT FACE) TO THE GRADIENTS
  m_faceTermComputers[currFaceIdx]->computeGradientFaceTerm(m_gradUpdates);

  // add the contribution to the gradients
  const CFuint currCellSide = (*m_currCellSide)[currFaceIdx];
  const CFuint nghbCellSide = currCellSide == LEFT ? RIGHT : LEFT;
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_nghbCellGrads[iSol])[iEq] += m_gradUpdates[nghbCellSide][iSol][iEq];
    }
  }

  // DIVIDE BY SOLUTION POINT JACOBIAN DETERMINANT
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDetNghbrCell[iSol];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_nghbCellGrads[iSol])[iEq] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::setup()
{
  CFAUTOTRACE;

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // maximum number of solution points and faces in a cell
  CFuint maxNbrSolPnts = 0;
  CFuint maxNbrFaces = 0;
  for (CFuint iElemType = 0; iElemType < sdLocalData.size(); ++iElemType)
  {
    const CFuint nbrSolPnts = sdLocalData[iElemType]->getNbrOfSolPnts();
    maxNbrSolPnts = maxNbrSolPnts > nbrSolPnts ? maxNbrSolPnts: nbrSolPnts;

    const CFuint nbrFaces = sdLocalData[iElemType]->getNbrCellFaces();
    maxNbrFaces = maxNbrFaces > nbrFaces ? maxNbrFaces : nbrFaces;
  }

  // get the number of equations and dimensionality
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim();

  // get cell builder
  m_cellBuilder = getMethodData().getCellBuilder();
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSide    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSide     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrients      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdx        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  // get the volume term computer
  m_volTermComputer = getMethodData().getVolTermComputer();

  // get face term computers and boundary face term computers
  m_faceTermComputers   .resize(maxNbrFaces);
  m_bndFaceTermComputers.resize(maxNbrFaces);
  m_faceTermComputers   [0] = getMethodData().getFaceTermComputer   ();
  m_bndFaceTermComputers[0] = getMethodData().getBndFaceTermComputer();
  CFuint faceCompIdx = 0;
  for (CFuint iFace = 1; iFace < maxNbrFaces; ++iFace, ++faceCompIdx)
  {
    m_faceTermComputers   [iFace] = getMethodData().getAdditionalFaceTermComputer   (faceCompIdx);
    m_bndFaceTermComputers[iFace] = getMethodData().getAdditionalBndFaceTermComputer(faceCompIdx);
  }

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

  // resize m_faceResUpdates and m_updateCoefFaceUpd
  m_faceResUpdates.resize(2);
  m_updateCoefFaceUpd.resize(2);

  // initialize to 10, to make sure that volume and face term computers' data is set initially
  m_currIter = 10;

  if (getMethodData().hasDiffTerm())
  {
    // get neighbouring cell builder
    m_nghbrCellBuilder = getMethodData().getSecondCellBuilder();
    m_isFaceOnBoundaryNghbrCell = m_nghbrCellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
    m_nghbrCellSideNghbrCell    = m_nghbrCellBuilder->getGeoBuilder()->getNeighbrCellSide ();
    m_currCellSideNghbrCell     = m_nghbrCellBuilder->getGeoBuilder()->getCurrentCellSide ();
    m_faceOrientsNghbrCell      = m_nghbrCellBuilder->getGeoBuilder()->getFaceOrient      ();
    m_faceBCIdxNghbrCell        = m_nghbrCellBuilder->getGeoBuilder()->getFaceBCIdx       ();

    // get the volume term compute for the neighbouring cell
    m_volTermComputerNghbrCell = getMethodData().getSecondVolTermComputer();

    // get face term computers and boundary face term computers for the neighbouring cells
    const CFuint maxNbrFacesM1 = maxNbrFaces - 1;
    m_faceTermComputersNghbrCell   .resize(maxNbrFacesM1);
    m_bndFaceTermComputersNghbrCell.resize(maxNbrFacesM1);
    for (CFuint iFace = 0; iFace < maxNbrFacesM1; ++iFace, ++faceCompIdx)
    {
      m_faceTermComputersNghbrCell   [iFace] = getMethodData().getAdditionalFaceTermComputer   (faceCompIdx);
      m_bndFaceTermComputersNghbrCell[iFace] = getMethodData().getAdditionalBndFaceTermComputer(faceCompIdx);
    }

    // resize m_faceNghbrStatesNghbrCell
    m_faceNghbrStatesNghbrCell.resize(maxNbrFacesM1);
    for (CFuint iFace = 0; iFace < maxNbrFacesM1; ++iFace)
    {
      m_faceNghbrStatesNghbrCell[iFace].resize(2);
    }

    // resize m_otherFaceLocalIdxs
    m_otherFaceLocalIdxs.resize(maxNbrFacesM1);

    // resize m_gradUpdates
    m_gradUpdates.resize(2);

    // resize m_solJacobDet and m_solJacobDetNghbrCell
    m_solJacobDet         .resize(maxNbrSolPnts);
    m_solJacobDetNghbrCell.resize(maxNbrSolPnts);

    // allocate memory for m_currCellGrads and m_nghbCellGrads
    m_currCellGrads.resize(maxNbrSolPnts);
    m_nghbCellGrads.resize(maxNbrSolPnts);
    for (CFuint iSol = 0; iSol < maxNbrSolPnts; ++iSol)
    {
      m_currCellGrads[iSol] = new vector<RealVector>(m_nbrEqs);
      m_nghbCellGrads[iSol] = new vector<RealVector>(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        (*m_currCellGrads[iSol])[iEq].resize(m_dim);
        (*m_nghbCellGrads[iSol])[iEq].resize(m_dim);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellSpectralFD::unsetup()
{
  CFAUTOTRACE;

  cf_assert(m_nghbCellGrads.size() == m_currCellGrads.size());
  for (CFuint iSol = 0; iSol < m_currCellGrads.size(); ++iSol)
  {
    deletePtr(m_currCellGrads[iSol]);
    deletePtr(m_nghbCellGrads[iSol]);
  }
  m_currCellGrads.resize(0);
  m_nghbCellGrads.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

vector< Common::SafePtr< BaseDataSocketSink > >
    RhsInGivenCellSpectralFD::needsSockets()
{
  vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_rhsCurrStatesSet);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_statesSetIdx);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
