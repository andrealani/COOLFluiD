#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/RhsInGivenCellCompactSpectralFD.hh"
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

MethodCommandProvider<RhsInGivenCellCompactSpectralFD, SpectralFDMethodData, SpectralFDModule>
    RhsInGivenCellCompactSpectralFDProvider("RhsInGivenCellCompact");

//////////////////////////////////////////////////////////////////////////////

RhsInGivenCellCompactSpectralFD::RhsInGivenCellCompactSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_updateCoeff("updateCoeff"),
  socket_statesSetIdx("statesSetIdx"),
  m_cellBuilder(CFNULL),
  m_volTermComputer(CFNULL),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_bcStateComputers(),
  m_cell(CFNULL),
  m_faces(CFNULL),
  m_cellStates(CFNULL),
  m_faceNghbrStates(),
  m_isFaceOnBoundary(CFNULL),
  m_nghbrCellSide(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_resUpdates(),
  m_diffResUpdates(),
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
  m_solJacobDet(),
  m_solPntsLocalCoords(CFNULL),
  m_otherFaceLocalIdxs()
{
}

//////////////////////////////////////////////////////////////////////////////

RhsInGivenCellCompactSpectralFD::~RhsInGivenCellCompactSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellCompactSpectralFD::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellCompactSpectralFD::execute()
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
    // get the neighbouring faces
    m_faces = m_cell->getNeighborGeos();

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

    // COMPUTE THE FACE TERMS
    setFaceTermDataAndComputeFaceTermsAndUpdateCoefs();

    // COMPUTE THE VOLUME TERM
    setVolumeTermDataAndComputeCellGradientsAndVolumeTerm();

    // MULTIPLY RHS WITH RESIDUAL FACTOR
    multplyRHSWithResFactor();
  }

  //release the GeometricEntity
  m_cellBuilder->releaseGE();
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellCompactSpectralFD::setVolumeAndFaceTermComputersData()
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

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get solution point local coordinates
  m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellCompactSpectralFD::resizeResUpdates()
{
  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // resize m_resUpdates
  const CFuint nbrRes = nbrSolPnts*m_nbrEqs;
  m_resUpdates    .resize(nbrRes);
  m_diffResUpdates.resize(nbrRes);
  m_faceResUpdates    [LEFT ].resize(nbrRes);
  m_faceResUpdates    [RIGHT].resize(nbrRes);

  // resize m_gradUpdates
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

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellCompactSpectralFD::clearResidual()
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

void RhsInGivenCellCompactSpectralFD::addUpdatesToResidual()
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

void RhsInGivenCellCompactSpectralFD::addUpdateToUpdateCoef()
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

void RhsInGivenCellCompactSpectralFD::setFaceTermDataAndComputeFaceTermsAndUpdateCoefs()
{
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    // set face neighbour states
    const CFuint nbrFaceNeighbours = (*m_faces)[iFace]->nbNeighborGeos();
    for (CFuint iSide = 0; iSide < nbrFaceNeighbours; ++iSide)
    {
      GeometricEntity* cell = (*m_faces)[iFace]->getNeighborGeo(iSide);
      m_faceNghbrStates[iFace][iSide] = cell->getStates();
    }

    if ((*m_isFaceOnBoundary)[iFace])
    {
      // set BC state computer
      m_bndFaceTermComputers[iFace]->setBcStateComputer((*m_bcStateComputers)[(*m_faceBCIdx)[iFace]]);

      // set the orientation of the face (should be set here, might be changed in the command that computes the face terms)
      m_bndFaceTermComputers[iFace]->setFaceOrientation(iFace);

      // set the face in the boundary face term computer
      m_bndFaceTermComputers[iFace]->setCurrentFace((*m_faces)[iFace]);

      // compute the face data in the boundary face term computer
      m_bndFaceTermComputers[iFace]->computeFaceData();

      // compute the neighbouring cell data in the boundary face term computer
      m_bndFaceTermComputers[iFace]->computeNeighbourCellData();

      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->reconstructFluxPntsStates(*m_cellStates);

      // compute the solution polynomial gradients in the face flux points
      m_bndFaceTermComputers[iFace]->reconstructFluxPntsSolPolyGrads(*m_cellStates);

      // compute the face term
      if (m_computeUpdateCoef)
      {
        m_bndFaceTermComputers[iFace]
          ->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_updateCoefUpd);
        m_bndFaceTermComputers[iFace]
          ->computeDiffFaceTermAndUpdateCoefContributions(m_diffResUpdates,m_updateCoefFaceUpd[0]);
        m_resUpdates += m_diffResUpdates;
        m_updateCoefUpd += m_updateCoefFaceUpd[0];
      }
      else
      {
        m_bndFaceTermComputers[iFace]->computeConvFaceTerm(m_resUpdates);
        m_bndFaceTermComputers[iFace]->computeDiffFaceTerm(m_diffResUpdates);
        m_resUpdates += m_diffResUpdates;
      }
    }
    else
    {
      // set the orientation of the face
      m_faceTermComputers[iFace]->setFaceOrientation((*m_faceOrients)[iFace]);

      // set the face in the face term computer
      m_faceTermComputers[iFace]->setCurrentFace((*m_faces)[iFace]);

      // compute the face data in the face term computer
      m_faceTermComputers[iFace]->computeFaceData();

      // compute the neighbouring cell data in the face term computer
      m_faceTermComputers[iFace]->computeNeighbourCellData();

      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->reconstructFluxPntsStates(m_faceNghbrStates[iFace]);

      // compute the solution polynomial gradients in the face flux points
      m_faceTermComputers[iFace]->reconstructFluxPntsSolPolyGrads(m_faceNghbrStates[iFace]);

      // get current cell side
      const CFuint currCellSide = (*m_currCellSide)[iFace];

      // compute the convective face term
      if (m_computeUpdateCoef)
      {
        m_faceTermComputers[iFace]
            ->computeConvFaceTermAndWaveSpeedUpdates(m_faceResUpdates,m_updateCoefFaceUpd);
        m_resUpdates = m_faceResUpdates[currCellSide];
        m_updateCoefUpd = m_updateCoefFaceUpd[currCellSide];
        m_faceTermComputers[iFace]
            ->computeDiffFaceTermAndUpdateCoefContributions(m_faceResUpdates,m_updateCoefFaceUpd);
        m_resUpdates += m_faceResUpdates[currCellSide];
        m_updateCoefUpd += m_updateCoefFaceUpd[currCellSide];
      }
      else
      {
        m_faceTermComputers[iFace]->computeConvFaceTerm(m_faceResUpdates);
        m_resUpdates = m_faceResUpdates[currCellSide];
        m_faceTermComputers[iFace]->computeDiffFaceTerm(m_faceResUpdates);
        m_resUpdates += m_faceResUpdates[currCellSide];
      }
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

void RhsInGivenCellCompactSpectralFD::setVolumeTermDataAndComputeCellGradientsAndVolumeTerm()
{
  // set the current cell and compute the cell data in the volume term computer
  m_volTermComputer->setCurrentCell(m_cell);
  m_volTermComputer->computeCellData();

  // reconstruct the solution in the flux points
  m_volTermComputer->reconstructStates(*m_cellStates);

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

  // compute the volume term
  m_volTermComputer->computeCellConvVolumeTerm(m_resUpdates);
  m_volTermComputer->computeCellDiffVolumeTerm(m_diffResUpdates);
  m_resUpdates += m_diffResUpdates;

  // add volume term to the current cell rhs
  addUpdatesToResidual();
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellCompactSpectralFD::multplyRHSWithResFactor()
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

void RhsInGivenCellCompactSpectralFD::setup()
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
  m_volTermComputer = getMethodData().getVolTermComputer().d_castTo<CompactVolTermComputer>();

  // get face term computers and boundary face term computers
  m_faceTermComputers   .resize(maxNbrFaces);
  m_bndFaceTermComputers.resize(maxNbrFaces);
  m_faceTermComputers   [0] = getMethodData().getFaceTermComputer   ().d_castTo<CompactFaceTermComputer   >();
  m_bndFaceTermComputers[0] = getMethodData().getBndFaceTermComputer().d_castTo<CompactBndFaceTermComputer>();
  CFuint faceCompIdx = 0;
  for (CFuint iFace = 1; iFace < maxNbrFaces; ++iFace, ++faceCompIdx)
  {
    m_faceTermComputers   [iFace] = getMethodData().
        getAdditionalFaceTermComputer   (faceCompIdx).d_castTo<CompactFaceTermComputer   >();
    m_bndFaceTermComputers[iFace] = getMethodData().
        getAdditionalBndFaceTermComputer(faceCompIdx).d_castTo<CompactBndFaceTermComputer>();
  }

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

  // resize variables
  m_faceResUpdates   .resize(2);
  m_updateCoefFaceUpd.resize(2);

  // initialize to 10, to make sure that volume and face term computers' data is set initially
  m_currIter = 10;

  // resize m_gradUpdates
  m_gradUpdates.resize(2);

  // resize m_solJacobDet and m_solJacobDetNghbrCell
  m_solJacobDet.resize(maxNbrSolPnts);

  // allocate memory for m_currCellGrads
  m_currCellGrads.resize(maxNbrSolPnts);
  for (CFuint iSol = 0; iSol < maxNbrSolPnts; ++iSol)
  {
    m_currCellGrads[iSol] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_currCellGrads[iSol])[iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RhsInGivenCellCompactSpectralFD::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSol = 0; iSol < m_currCellGrads.size(); ++iSol)
  {
    deletePtr(m_currCellGrads[iSol]);
  }
  m_currCellGrads.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    RhsInGivenCellCompactSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_rhsCurrStatesSet);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_statesSetIdx);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
