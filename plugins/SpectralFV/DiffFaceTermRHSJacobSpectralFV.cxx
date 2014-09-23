#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFV/DiffFaceTermRHSJacobSpectralFV.hh"
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

MethodCommandProvider<DiffFaceTermRHSJacobSpectralFV, SpectralFVMethodData, SpectralFVModule> DiffFaceTermRHSJacobSpectralFVProvider("DiffFaceTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

DiffFaceTermRHSJacobSpectralFV::DiffFaceTermRHSJacobSpectralFV(const std::string& name) :
  DiffFaceTermRHSSpectralFV(name),
  m_cellBuilders(),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_accSC(CFNULL),
  m_volTermComputers(),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_faces(),
  m_faceNghbrStates(),
  m_faceNghbrGrads(),
  m_invVolFracCVs(),
  m_svFaceCVConn(),
  m_pertResUpdates(),
  m_derivResUpdates(),
  m_gradUpdates(),
  m_pertGrads(),
  m_cellGradsMinusFaceTerm(),
  m_cellGradsMinusOtherFaceTerm(),
  m_unpertCellDiffRes(),
  m_pertCellDiffRes(),
  m_derivCellDiffRes(),
  m_cvInvVols(),
  m_otherFaceLocalIdxs(),
  m_isFaceOnBoundary(),
  m_nghbrCellSide(),
  m_currCellSide(),
  m_faceOrients(),
  m_faceBCIdx(),
  m_bcStateComputers(CFNULL)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

DiffFaceTermRHSJacobSpectralFV::~DiffFaceTermRHSJacobSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::configure ( Config::ConfigArgs& args )
{
  DiffFaceTermRHSSpectralFV::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::execute()
{
  CFTRACEBEGIN;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

  // set the data needed to compute the face terms;
  setFaceTermData();

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

        // get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current boundary faces
        setOtherFacesLocalIdxs();

        // get the neigbouring states of the other faces
        setFaceNeighbourStates();

        // set the gradients
        setGradients();

        // get the neigbouring gradients of the other faces
        setFaceNeighbourGradients();

        // set the current face and compute the face data in the face term computer
        m_faceTermComputer->setCurrentFace(m_face);
        m_faceTermComputer->computeFaceData();
        m_faceTermComputer->computeNeighbourCellData();

        // set the neighbouring cells data
        setCellsData();

        // reconstruct the states
        m_faceTermComputer->reconstructFluxPntsStates(*m_states[LEFT],*m_states[RIGHT]);
        m_faceTermComputer->reconstructFaceAvgState  (*m_states[LEFT],*m_states[RIGHT]);

        // reconstruct the gradients
        m_faceTermComputer->reconstructFluxPntsGradients(m_grads[LEFT],m_grads[RIGHT],
                                                        m_grads[LEFT].size(),m_grads[RIGHT].size());

        // compute the face terms and the update coefficient contributions
        m_faceTermComputer
                    ->computeDiffFaceTermAndUpdateCoefContributions(m_resUpdates,
                                                                    m_updateCoefContr[LEFT ],
                                                                    m_updateCoefContr[RIGHT]);

        // update the rhs
        updateRHS();

        // add update coefficient contributions
        addUpdateCoeffContributions();

        // compute auxiliary term for the perturbed gradient reconstruction from current face
        computeCellGradsMinusFaceTerm();

        // compute the unperturbed cell diffusive residuals
        computeUnpertCellDiffResiduals();
      }

      // compute the diffusive face term contribution to the jacobian
      if (lParUpdatable && rParUpdatable)
      {
        // compute auxiliary terms for the perturbed gradient reconstruction from other faces
        computeCellGradsMinusOtherFaceTerms(LEFT );
        computeCellGradsMinusOtherFaceTerms(RIGHT);

        computeBothJacobsDiffFaceTerm();
      }
      else if (lParUpdatable)
      {
        // compute auxiliary terms for the perturbed gradient reconstruction from other faces
        computeCellGradsMinusOtherFaceTerms(RIGHT);

        computeOneJacobDiffFaceTerm(LEFT );
      }
      else if (rParUpdatable)
      {
        // compute auxiliary terms for the perturbed gradient reconstruction from other faces
        computeCellGradsMinusOtherFaceTerms(LEFT );

        computeOneJacobDiffFaceTerm(RIGHT);
      }

      // release the cells
      if (lParUpdatable || rParUpdatable)
      {
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
      }

      // release the face
      m_faceBuilder->releaseGE();
    }
  }

// CF_DEBUG_EXIT;
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::setFaceTermData()
{
  DiffFaceTermRHSSpectralFV::setFaceTermData();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get the volume fractions of the CVs
  m_invVolFracCVs = svLocalData[0]->getInvVolFracCV();

  // get the boundary SV face-CV connectivity
  m_svFaceCVConn = svLocalData[0]->getExtSVFaceCVConn();

  // set the data in the face term computers
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iFace = 0; iFace < m_faceTermComputers[iSide].size(); ++iFace)
    {
      m_faceTermComputers[iSide][iFace]->setFaceTermData();
    }
  }

  // set the data in the boundary face term computers
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iFace = 0; iFace < m_bndFaceTermComputers[iSide].size(); ++iFace)
    {
      m_bndFaceTermComputers[iSide][iFace]->setFaceTermData();
    }
  }

  // set the data in the volume term computers
  m_volTermComputers[LEFT ]->setVolumeTermData(0);
  m_volTermComputers[RIGHT]->setVolumeTermData(0);
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::setOtherFacesLocalIdxs()
{
  // get face ID of current boundary face
  const CFuint currFaceID = m_face->getID();

  // loop over the faces of the cells
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrFaces = m_cells[iSide]->nbNeighborGeos();
    cf_assert(m_faces[iSide]->size() == nbrFaces);
    CFuint iFace = 0;
    for (CFuint faceIdx = 0; faceIdx < nbrFaces; ++faceIdx)
    {
      if ((*m_faces[iSide])[faceIdx]->getID() != currFaceID)
      {
        m_otherFaceLocalIdxs[iSide][iFace] = faceIdx;
        ++iFace;
      }
    }
    cf_assert(iFace == m_otherFaceLocalIdxs[iSide].size());
    cf_assert(iFace == nbrFaces-1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::setFaceNeighbourStates()
{
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbour states of other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      // get neigbouring states
      const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundary[iSide])[faceIdx] ? 1 : 2;
      for (CFuint iSide2 = 0; iSide2 < nbrFaceNeighbours; ++iSide2)
      {
        GeometricEntity* cell = (*m_faces[iSide])[faceIdx]->getNeighborGeo(iSide2);
        m_faceNghbrStates[iSide][iFace][iSide2] = cell->getStates();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::setFaceNeighbourGradients()
{
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbour states of other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      // get neigbouring states
      const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundary[iSide])[faceIdx] ? 1 : 2;
      for (CFuint iSide2 = 0; iSide2 < nbrFaceNeighbours; ++iSide2)
      {
        // get number of states
        const CFuint nbrStates = m_faceNghbrStates[iSide][iFace][iSide2]->size();

        // resize m_faceNghbrGrads[iSide][iFace][iSide2]
        m_faceNghbrGrads[iSide][iFace][iSide2].resize(nbrStates);

        // set the gradients
        for (CFuint iState = 0; iState < nbrStates; ++iState)
        {
          const CFuint stateID = (*m_faceNghbrStates[iSide][iFace][iSide2])[iState]->getLocalID();
          m_faceNghbrGrads[iSide][iFace][iSide2][iState] = &gradients[stateID];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::setCellsData()
{
  // set cells data
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // cell volume
    const CFreal cellInvVol = 1.0/m_cells[iSide]->computeVolume();

    // cv volumes
    const CFuint nbrCVs = m_invVolFracCVs->size();
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
    {
      m_cvInvVols[iSide][iCV] = cellInvVol*(*m_invVolFracCVs)[iCV];
    }

    // set the neighbouring cell and compute neighbouring cell data in the volume term computer
    m_volTermComputers[iSide]->setCurrentCell(m_cells[iSide]);
    m_volTermComputers[iSide]->computeCellData();
    m_volTermComputers[iSide]->reconstructStates(*m_states[iSide]);
    m_volTermComputers[iSide]->reconstructGradients(m_grads[iSide],m_grads[iSide].size());

    // set orientation or boundary condition state computers of the other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      // set orientation or boundary condition
      if ((*m_isFaceOnBoundary[iSide])[faceIdx])
      {
        // set BC state computer
        m_bndFaceTermComputers[iSide][iFace]
            ->setBcStateComputer((*m_bcStateComputers)[(*m_faceBCIdx[iSide])[faceIdx]]);

        // set the orientation of the face
        m_bndFaceTermComputers[iSide][iFace]->setFaceOrientation(faceIdx);

        // set the face in the boundary face term computer
        m_bndFaceTermComputers[iSide][iFace]->setCurrentFace((*m_faces[iSide])[faceIdx]);

        // compute the face data in the boundary face term computer
        m_bndFaceTermComputers[iSide][iFace]->computeFaceData();

        // compute the neighbouring cell data in the boundary face term computer
        m_bndFaceTermComputers[iSide][iFace]->computeNeighbourCellData();

        // reconstruct the states in the face flux points
        m_bndFaceTermComputers[iSide][iFace]->reconstructFluxPntsStates(*m_states[iSide]);

        // reconstruct the gradients in the face flux points
        m_bndFaceTermComputers[iSide][iFace]->
            reconstructFluxPntsGradients(m_grads[iSide],m_grads[iSide].size());
      }
      else
      {
        // set the orientation of the face
        m_faceTermComputers[iSide][iFace]->setFaceOrientation((*m_faceOrients[iSide])[faceIdx]);

        // set the face in the face term computer
        m_faceTermComputers[iSide][iFace]->setCurrentFace((*m_faces[iSide])[faceIdx]);

        // compute the face data in the face term computer
        m_faceTermComputers[iSide][iFace]->computeFaceData();

        // compute the neighbouring cell data in the face term computer
        m_faceTermComputers[iSide][iFace]->computeNeighbourCellData();

        // reconstruct the states in the face flux points
        m_faceTermComputers[iSide][iFace]->
            reconstructFluxPntsStates(*m_faceNghbrStates[iSide][iFace][LEFT ],
                                      *m_faceNghbrStates[iSide][iFace][RIGHT]);

        // reconstruct the gradients in the face flux points
        m_faceTermComputers[iSide][iFace]->
            reconstructFluxPntsGradients(m_faceNghbrGrads[iSide][iFace][LEFT ],
                                         m_faceNghbrGrads[iSide][iFace][RIGHT],
                                         m_faceNghbrGrads[iSide][iFace][LEFT ].size(),
                                         m_faceNghbrGrads[iSide][iFace][RIGHT].size());
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computeCellGradsMinusFaceTerm()
{
  // copy cell gradients
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrCVs = m_grads[iSide].size();
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_cellGradsMinusFaceTerm[iSide][iCV][iEq] = (*m_grads[iSide][iCV])[iEq];
      }
    }
  }

  // compute current face contribution to the gradients
  m_faceTermComputer->computeGradientFaceTerm(m_gradUpdates);

  // subtract the face term
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get internal face neighbour CV index
    const CFuint cvIdxL  = (*m_cvCVConnSVFace)[m_orient][iSubFace][LEFT ];
    const CFuint cvIdxR  = (*m_cvCVConnSVFace)[m_orient][iSubFace][RIGHT];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellGradsMinusFaceTerm[LEFT ][cvIdxL][iEq] -=
          m_gradUpdates[iSubFace][iEq]*m_cvInvVols[LEFT ][cvIdxL];
      m_cellGradsMinusFaceTerm[RIGHT][cvIdxR][iEq] +=
          m_gradUpdates[iSubFace][iEq]*m_cvInvVols[RIGHT][cvIdxR];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computeCellGradsMinusOtherFaceTerms(const CFuint side)
{
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];
    if(!(*m_isFaceOnBoundary[side])[faceIdx])
    {
      // copy cell gradients
      const CFuint nbrCVs = m_grads[side].size();
      for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_cellGradsMinusOtherFaceTerm[side][iFace][iCV][iEq] =
              (*m_grads[side][iCV])[iEq];
        }
      }

      // compute the internal face contribution to the gradients
      m_faceTermComputers[side][iFace]->computeGradientFaceTerm(m_gradUpdates);

      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide[side])[faceIdx];

      // orientation of this face
      const CFuint faceOrient = (*m_faceOrients[side])[faceIdx];

      // subtract the gradient face term
      const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(cellSide));
      const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[faceOrient].size();
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
      {
        // get internal face neighbour CV index
        const CFuint cvIdx = (*m_cvCVConnSVFace)[faceOrient][iSubFace][cellSide];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_cellGradsMinusOtherFaceTerm[side][iFace][cvIdx][iEq] -=
              sideFactor*m_gradUpdates[iSubFace][iEq]*m_cvInvVols[side][cvIdx];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computeUnpertCellDiffResiduals()
{
  // compute the volume terms
  m_volTermComputers[LEFT ]->computeCellDiffVolumeTerm(m_unpertCellDiffRes[LEFT ]);
  m_volTermComputers[RIGHT]->computeCellDiffVolumeTerm(m_unpertCellDiffRes[RIGHT]);

  // add current face diffusive fluxes (m_resUpdates is set outside this function)
  CFuint iRes = 0;
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get first IDs of residuals in these CVs
    CFuint resIDL = m_nbrEqs*(*m_cvCVConnSVFace)[m_orient][iSubFace][LEFT ];
    CFuint resIDR = m_nbrEqs*(*m_cvCVConnSVFace)[m_orient][iSubFace][RIGHT];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++iRes, ++resIDL, ++resIDR)
    {
      m_unpertCellDiffRes[LEFT ][resIDL] += m_resUpdates[iRes];
      m_unpertCellDiffRes[RIGHT][resIDR] -= m_resUpdates[iRes];
    }
  }

  // add other face diffusive fluxes
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get number of other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      if ((*m_isFaceOnBoundary[iSide])[faceIdx])
      {
        // compute the boundary face contribution to the diffusive residuals
        // using m_pertResUpdates because the values stored in m_resUpdates should be preserved
        m_bndFaceTermComputers[iSide][iFace]->computeDiffFaceTerm(m_pertResUpdates);

        // add the contribution to the diffusive residuals
        CFuint iRes = 0;
        const CFuint nbrSubFaces = (*m_svFaceCVConn)[faceIdx].size();
        for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
        {
          // get boundary face neighbour CV index
          CFuint resID = m_nbrEqs*(*m_svFaceCVConn)[faceIdx][iSubFace];
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID, ++iRes)
          {
            m_unpertCellDiffRes[iSide][resID] += m_pertResUpdates[iRes];
          }
        }
      }
      else
      {
        // compute the internal face contribution to the diffusive residuals
        // using m_pertResUpdates because the values stored in m_resUpdates should be preserved
        m_faceTermComputers[iSide][iFace]->computeDiffFaceTerm(m_pertResUpdates);

        // cell side with respect to this face
        const CFuint cellSide = (*m_currCellSide[iSide])[faceIdx];

        // orientation of this face
        const CFuint faceOrient = (*m_faceOrients[iSide])[faceIdx];

        // add the contribution to the diffusive residuals
        const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(cellSide));
        CFuint iRes = 0;
        const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[faceOrient].size();
        for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
        {
          // get internal face neighbour CV index
          CFuint resID = m_nbrEqs*(*m_cvCVConnSVFace)[faceOrient][iSubFace][cellSide];
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID, ++iRes)
          {
            m_unpertCellDiffRes[iSide][resID] += sideFactor*m_pertResUpdates[iRes];
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::
    backupAndReconstructOtherFacesAndCellPhysVars(const CFuint side, const CFuint iVar)
{
  // face term computers
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];

    // backup and reconstruct physical variable
    if ((*m_isFaceOnBoundary[side])[faceIdx])
    {
      m_bndFaceTermComputers[side][iFace]->backupAndReconstructPhysVar(iVar,*m_states[side]);
    }
    else
    {
      m_faceTermComputers[side][iFace]->
          backupAndReconstructPhysVar((*m_currCellSide[side])[faceIdx],iVar,*m_states[side]);
    }
  }

  // volume term computer
  m_volTermComputers[side]->backupAndReconstructPhysVar(iVar,*m_states[side]);
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::
    restoreOtherFacesAndCellPhysVars(const CFuint side, const CFuint iVar)
{
  // face term computers
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];

    // restore physical variable
    if ((*m_isFaceOnBoundary[side])[faceIdx])
    {
      m_bndFaceTermComputers[side][iFace]->restorePhysVar(iVar);
    }
    else
    {
      m_faceTermComputers[side][iFace]->restorePhysVar((*m_currCellSide[side])[faceIdx],iVar);
    }
  }

  // volume term computer
  m_volTermComputers[side]->restorePhysVar(iVar);
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computePerturbedGradients(const CFuint side)
{
  // compute the internal contribution to the gradients
  m_volTermComputers[side]->computeGradientVolumeTerm(m_gradUpdates);
  const CFuint nbrCVs = m_grads[side].size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[side][iCV])[iEq] = m_gradUpdates[iCV][iEq];
    }
  }

  // compute other face contributions to the gradients
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];

    if ((*m_isFaceOnBoundary[side])[faceIdx])
    {
      // compute the boundary face contribution to the gradients
      m_bndFaceTermComputers[side][iFace]->computeGradientFaceTerm(m_gradUpdates);

      // add the contribution to the gradients
      const CFuint nbrSubFaces = (*m_svFaceCVConn)[faceIdx].size();
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
      {
        // get boundary face neighbour CV index
        const CFuint cvIdx  = (*m_svFaceCVConn)[faceIdx][iSubFace];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_pertGrads[side][cvIdx])[iEq] += m_gradUpdates[iSubFace][iEq];
        }
      }
    }
    else
    {
      // compute the internal face contribution to the gradients
      m_faceTermComputers[side][iFace]->computeGradientFaceTerm(m_gradUpdates);

      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide[side])[faceIdx];

      // orientation of this face
      const CFuint faceOrient = (*m_faceOrients[side])[faceIdx];

      // add the contribution to the gradients
      const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(cellSide));
      const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[faceOrient].size();
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
      {
        // get internal face neighbour CV index
        const CFuint cvIdx = (*m_cvCVConnSVFace)[faceOrient][iSubFace][cellSide];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_pertGrads[side][cvIdx])[iEq] += sideFactor*m_gradUpdates[iSubFace][iEq];
        }
      }
    }
  }

  // compute current face contribution to the gradients
  m_faceTermComputer->computeGradientFaceTerm(m_gradUpdates);

  // add the contribution to the gradients
  const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(side));
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
     // get internal face neighbour CV index
    const CFuint cvIdx  = (*m_cvCVConnSVFace)[m_orient][iSubFace][side];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[side][cvIdx])[iEq] += sideFactor*m_gradUpdates[iSubFace][iEq];
    }
  }

  // divide by CV volumes
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[side][iCV])[iEq] *= m_cvInvVols[side][iCV];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computePertGradsFromFaceTerm(const CFuint side)
{
  // copy m_cellGradsMinusFaceTerm into pertGrads
  const CFuint nbrCVs = m_grads[side].size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[side][iCV])[iEq] = m_cellGradsMinusFaceTerm[side][iCV][iEq];
    }
  }

  // current face contribution to the gradients should have been computed in computePerturbedGradients()!

  // subtract the face term
  const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(side));
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get internal face neighbour CV index
    const CFuint cvIdx  = (*m_cvCVConnSVFace)[m_orient][iSubFace][side];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[side][cvIdx])[iEq] +=
          + sideFactor*m_cvInvVols[side][cvIdx]*m_gradUpdates[iSubFace][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computePertGradsFromOtherFaceTerm(const CFuint side, const CFuint iFace)
{
  // copy m_cellGradsMinusOtherFaceTerm into pertGrads
  const CFuint nbrCVs = m_grads[side].size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[side][iCV])[iEq] = m_cellGradsMinusOtherFaceTerm[side][iFace][iCV][iEq];
    }
  }

  // get local face index
  const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];
  cf_assert(!(*m_isFaceOnBoundary[side])[faceIdx]);

  // compute the internal face contribution to the gradients
  m_faceTermComputers[side][iFace]->computeGradientFaceTerm(m_gradUpdates);

  // cell side with respect to this face
  const CFuint cellSide = (*m_currCellSide[side])[faceIdx];

  // orientation of this face
  const CFuint faceOrient = (*m_faceOrients[side])[faceIdx];

  // subtract the gradient face term
  const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(cellSide));
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[faceOrient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get internal face neighbour CV index
    const CFuint cvIdx = (*m_cvCVConnSVFace)[faceOrient][iSubFace][cellSide];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[side][cvIdx])[iEq] +=
          sideFactor*m_gradUpdates[iSubFace][iEq]*m_cvInvVols[side][cvIdx];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::reconstructOtherFacesAndCellGradients(const CFuint side)
{
    // face term computers
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];

    // backup and reconstruct physical variable
    if ((*m_isFaceOnBoundary[side])[faceIdx])
    {
      m_bndFaceTermComputers[side][iFace]->
          reconstructFluxPntsGradients(m_pertGrads[side],m_grads[side].size());
    }
    else
    {
      m_faceTermComputers[side][iFace]->
          reconstructFluxPntsGradients((*m_currCellSide[side])[faceIdx],
                                       m_pertGrads[side],m_grads[side].size());
    }
  }

  // volume term computer
  m_volTermComputers[side]->reconstructGradients(m_pertGrads[side],m_grads[side].size());
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computePertCellDiffResiduals(const CFuint side)
{
  // compute the volume terms
  m_volTermComputers[side]->computeCellDiffVolumeTerm(m_pertCellDiffRes);

  // add current face diffusive fluxes (m_pertResUpdates is set outside this function)
  const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(side));
  CFuint iRes = 0;
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get first IDs of residuals in these CVs
    CFuint resID = m_nbrEqs*(*m_cvCVConnSVFace)[m_orient][iSubFace][side];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++iRes, ++resID)
    {
      m_pertCellDiffRes[resID] += sideFactor*m_pertResUpdates[iRes];
    }
  }

  // add other face diffusive fluxes
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];

    if ((*m_isFaceOnBoundary[side])[faceIdx])
    {
      // compute the boundary face contribution to the diffusive residuals
      m_bndFaceTermComputers[side][iFace]->computeDiffFaceTerm(m_pertResUpdates);

      // add the contribution to the diffusive residuals
      CFuint iRes = 0;
      const CFuint nbrSubFaces = (*m_svFaceCVConn)[faceIdx].size();
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
      {
        // get boundary face neighbour CV index
        CFuint resID = m_nbrEqs*(*m_svFaceCVConn)[faceIdx][iSubFace];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID, ++iRes)
        {
          m_pertCellDiffRes[resID] += m_pertResUpdates[iRes];
        }
      }
    }
    else
    {
      // compute the internal face contribution to the diffusive residuals
      m_faceTermComputers[side][iFace]->computeDiffFaceTerm(m_pertResUpdates);

      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide[side])[faceIdx];

      // orientation of this face
      const CFuint faceOrient = (*m_faceOrients[side])[faceIdx];

      // add the contribution to the diffusive residuals
      const CFreal sideFactor = pow(-1.0,static_cast<CFreal>(cellSide));
      CFuint iRes = 0;
      const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[faceOrient].size();
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
      {
        // get internal face neighbour CV index
        CFuint resID = m_nbrEqs*(*m_cvCVConnSVFace)[faceOrient][iSubFace][cellSide];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID, ++iRes)
        {
          m_pertCellDiffRes[resID] += sideFactor*m_pertResUpdates[iRes];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computeBothJacobsDiffFaceTerm()
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
    // variable for the other side
    const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;

    // perturbed side factor
    const CFreal pertSideFac = resFactor*pow(-1.0,static_cast<CFreal>(iSide));

    // loop over the states/CVs to perturb the states
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

        // backup and reconstruct physical variable in the face flux points
        m_faceTermComputer->backupAndReconstructPhysVar(iSide,iEqPert,*m_states[iSide]);

        // backup and reconstruct physical variable in the other face flux points
        // and in the cell flux points
        backupAndReconstructOtherFacesAndCellPhysVars(iSide,iEqPert);

        // compute the perturbed gradients in the current cell
        computePerturbedGradients(iSide);

        // compute the perturbed gradients in the other cell
        computePertGradsFromFaceTerm(iOtherSide);

        // reconstruct the face gradients
        m_faceTermComputer->reconstructFluxPntsGradients(m_pertGrads[LEFT]   ,m_pertGrads[RIGHT],
                                                         m_grads[LEFT].size(),m_grads[RIGHT].size());

        // compute the perturbed face term
        m_faceTermComputer->computeDiffFaceTerm(m_pertResUpdates);

        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

        // multiply residual update derivatives with residual factor and +1.0 or -1.0 depending on iSide
        m_derivResUpdates *= pertSideFac;

        // add the derivative of the residual updates to the accumulator
        const CFuint pertSideTerm = iSide*nbrLeftCVs;
        CFuint resUpdIdx = 0;
        for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx +=m_nbrEqs)
        {
          CFuint cvIdx = (*m_cvCVConnSVFace)[m_orient][iSubFace][iSide];
          acc.addValues(cvIdx+pertSideTerm,iCVPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
        }

        // reconstruct the gradients for the other cell and the other cell neighbouring faces
        reconstructOtherFacesAndCellGradients(iOtherSide);

        // compute the perturbed diffusive residual in the other cell
        computePertCellDiffResiduals(iOtherSide);

        // compute the finite difference derivative of the other cell the diffusive residual
        m_numJacob->computeDerivative(m_pertCellDiffRes,m_unpertCellDiffRes[iOtherSide],m_derivCellDiffRes);

        // multiply residual update derivatives with residual factor
        m_derivCellDiffRes *= resFactor;

        // add the derivative of the residual updates to the accumulator
        const CFuint otherSideTerm = iOtherSide*nbrLeftCVs;
        resUpdIdx = 0;
        for (CFuint iCV = 0; iCV < nbrCVs; ++iCV, resUpdIdx += m_nbrEqs)
        {
          acc.addValues(iCV+otherSideTerm,iCVPert+pertSideTerm,iEqPert,&m_derivCellDiffRes[resUpdIdx]);
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[iEqPert]);

        // restore physical variable in the flux points (face term computer)
        m_faceTermComputer->restorePhysVar(iSide,iEqPert);

        // restore physical variable in the other face flux points
        // and in the cell flux points
        restoreOtherFacesAndCellPhysVars(iSide,iEqPert);
      }
    }
  }
//   acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();

  // dereference single cell accumulator
  BlockAccumulator& accSC = *m_accSC;

  // compute the contribution of neighbours of the face neighbour cells
  // loop over the two sides with respect to the face
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // variable for the other side
    const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;

    // number of CVs in current face neighbour cell
    const CFuint nbrNeighbrCellCVs = m_states[iOtherSide]->size();

    // perturbed side factor
    const CFreal pertSideFac = resFactor*pow(-1.0,static_cast<CFreal>(iOtherSide));

    // reconstruct the other face gradients using the unperturbed gradients
    m_faceTermComputer->
        reconstructFluxPntsGradients(iOtherSide,m_grads[iOtherSide],m_grads[iOtherSide].size());

    // get number of other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      // if inner face, compute the additional jacobian contributions
      if (!(*m_isFaceOnBoundary[iSide])[faceIdx])
      {
        // get contributing cell side with respect to face
        const CFuint contrCellSide = (*m_nghbrCellSide[iSide])[faceIdx];

        // get contributing cell states
        vector< State* >* contrCellStates = m_faceNghbrStates[iSide][iFace][contrCellSide];
        const CFuint contrCellNbrStates = contrCellStates->size();

        // set block row and column indices
        // neighbour cell
        for (CFuint iCV = 0; iCV < nbrNeighbrCellCVs; ++iCV)
        {
          accSC.setRowIndex(iCV,(*m_states[iOtherSide])[iCV]->getLocalID());
        }
        // contributing cell
        for (CFuint iCV = 0; iCV < contrCellNbrStates; ++iCV)
        {
          accSC.setColIndex(iCV,(*contrCellStates)[iCV]->getLocalID());
        }

        // get contributing face term computer
        SafePtr< BaseFaceTermComputer >
            contrFaceTermComputer = m_faceTermComputers[iSide][iFace];

        // loop over the states/CVs in the contributing cell to perturb the states
        for (CFuint iCVPert = 0; iCVPert < contrCellNbrStates; ++iCVPert)
        {
          // dereference state
          State& pertState = *(*contrCellStates)[iCVPert];

          // loop over the variables in the state
          for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
          {
            // perturb physical variable in state
            m_numJacob->perturb(iEqPert,pertState[iEqPert]);

            // backup and reconstruct physical variable in the contributing face flux points
            contrFaceTermComputer->
                backupAndReconstructPhysVar(contrCellSide,iEqPert,*contrCellStates);

            // compute the perturbed gradients in the cell on side iSide
            computePertGradsFromOtherFaceTerm(iSide,iFace);

            // reconstruct the current face gradients
            m_faceTermComputer->
                reconstructFluxPntsGradients(iSide,m_pertGrads[iSide],m_grads[iSide].size());

            // compute the perturbed face term
            m_faceTermComputer->computeDiffFaceTerm(m_pertResUpdates);

            // compute the finite difference derivative of the face term
            m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

            // multiply residual update derivatives with residual factor and +1.0 or -1.0 depending on iSide
            m_derivResUpdates *= pertSideFac;

            // add the derivative of the residual updates to the accumulator
            CFuint resUpdIdx = 0;
            for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx += m_nbrEqs)
            {
              CFuint cvIdx = (*m_cvCVConnSVFace)[m_orient][iSubFace][iOtherSide];
              accSC.addValues(cvIdx,iCVPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
            }

            // restore physical variable in state
            m_numJacob->restore(pertState[iEqPert]);

            // restore physical variable in the flux points (contributing face term computer)
            contrFaceTermComputer->restorePhysVar(contrCellSide,iEqPert);
          }
        }
//   accSC.printToScreen();

        if (getMethodData().doComputeJacobian())
        {
          // add the values to the jacobian matrix
          m_lss->getMatrix()->addValues(accSC);
        }

        // reset to zero the entries in the block accumulator
        accSC.reset();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::computeOneJacobDiffFaceTerm(const CFuint side)
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
    // variable for the other side
    const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;

    // perturbed side factor
    const CFreal pertSideFac = resFactor*pow(-1.0,static_cast<CFreal>(iSide));

    // loop over the states/CVs to perturb the states
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

        // backup and reconstruct physical variable in the face flux points
        m_faceTermComputer->backupAndReconstructPhysVar(iSide,iEqPert,*m_states[iSide]);

        // backup and reconstruct physical variable in the other face flux points
        // and in the cell flux points
        backupAndReconstructOtherFacesAndCellPhysVars(iSide,iEqPert);

        // compute the perturbed gradients in the current cell
        computePerturbedGradients(iSide);

        // compute the perturbed gradients in the other cell
        computePertGradsFromFaceTerm(iOtherSide);

        // reconstruct the face gradients
        m_faceTermComputer->reconstructFluxPntsGradients(m_pertGrads[LEFT]   ,m_pertGrads[RIGHT],
                                                         m_grads[LEFT].size(),m_grads[RIGHT].size());

        // compute the perturbed face term
        m_faceTermComputer->computeDiffFaceTerm(m_pertResUpdates);

        if (iSide == side)
        {
          // compute the finite difference derivative of the face term
          m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

          // multiply residual update derivatives with residual factor and +1.0 or -1.0 depending on iSide
          m_derivResUpdates *= pertSideFac;

          // add the derivative of the residual updates to the accumulator
          const CFuint pertSideTerm = iSide*nbrLeftCVs;
          CFuint resUpdIdx = 0;
          for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx +=m_nbrEqs)
          {
            CFuint cvIdx = (*m_cvCVConnSVFace)[m_orient][iSubFace][iSide];
            acc.addValues(cvIdx+pertSideTerm,iCVPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
          }
        }
        else
        {
          // reconstruct the gradients for the other cell and the other cell neighbouring faces
          reconstructOtherFacesAndCellGradients(iOtherSide);

          // compute the perturbed diffusive residual in the other cell
          computePertCellDiffResiduals(iOtherSide);

          // compute the finite difference derivative of the other cell the diffusive residual
          m_numJacob->computeDerivative(m_pertCellDiffRes,m_unpertCellDiffRes[iOtherSide],m_derivCellDiffRes);

          // multiply residual update derivatives with residual factor
          m_derivCellDiffRes *= resFactor;

          // add the derivative of the residual updates to the accumulator
          const CFuint pertSideTerm  = iSide*nbrLeftCVs;
          const CFuint otherSideTerm = iOtherSide*nbrLeftCVs;
          CFuint resUpdIdx = 0;
          for (CFuint iCV = 0; iCV < nbrCVs; ++iCV, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iCV+otherSideTerm,iCVPert+pertSideTerm,iEqPert,&m_derivCellDiffRes[resUpdIdx]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[iEqPert]);

        // restore physical variable in the flux points (face term computer)
        m_faceTermComputer->restorePhysVar(iSide,iEqPert);

        // restore physical variable in the other face flux points
        // and in the cell flux points
        restoreOtherFacesAndCellPhysVars(iSide,iEqPert);
      }
    }
  }
//   acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();

  // dereference single cell accumulator
  BlockAccumulator& accSC = *m_accSC;

  // compute the contribution of neighbours of the face neighbour cells
  // variable for the other side
  const CFuint otherSide = side == LEFT ? RIGHT : LEFT;

  // number of CVs in current face neighbour cell
  const CFuint nbrNeighbrCellCVs = m_states[side]->size();

  // perturbed side factor
  const CFreal pertSideFac = resFactor*pow(-1.0,static_cast<CFreal>(side));

  // reconstruct the other face gradients using the unperturbed gradients
  m_faceTermComputer->
      reconstructFluxPntsGradients(side,m_grads[side],m_grads[side].size());

  // get number of other faces
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[otherSide].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[otherSide][iFace];

    // if inner face, compute the additional jacobian contributions
    if (!(*m_isFaceOnBoundary[otherSide])[faceIdx])
    {
      // get contributing cell side with respect to face
      const CFuint contrCellSide = (*m_nghbrCellSide[otherSide])[faceIdx];

      // get contributing cell states
      vector< State* >* contrCellStates = m_faceNghbrStates[otherSide][iFace][contrCellSide];
      const CFuint contrCellNbrStates = contrCellStates->size();

      // set block row and column indices
      // neighbour cell
      for (CFuint iCV = 0; iCV < nbrNeighbrCellCVs; ++iCV)
      {
        accSC.setRowIndex(iCV,(*m_states[side])[iCV]->getLocalID());
      }
      // contributing cell
      for (CFuint iCV = 0; iCV < contrCellNbrStates; ++iCV)
      {
        accSC.setColIndex(iCV,(*contrCellStates)[iCV]->getLocalID());
      }

      // get contributing face term computer
      SafePtr< BaseFaceTermComputer >
          contrFaceTermComputer = m_faceTermComputers[otherSide][iFace];

      // loop over the states/CVs in the contributing cell to perturb the states
      for (CFuint iCVPert = 0; iCVPert < contrCellNbrStates; ++iCVPert)
      {
        // dereference state
        State& pertState = *(*contrCellStates)[iCVPert];

        // loop over the variables in the state
        for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
        {
          // perturb physical variable in state
          m_numJacob->perturb(iEqPert,pertState[iEqPert]);

          // backup and reconstruct physical variable in the contributing face flux points
          contrFaceTermComputer->
              backupAndReconstructPhysVar(contrCellSide,iEqPert,*contrCellStates);

          // compute the perturbed gradients in the cell on side iSide
          computePertGradsFromOtherFaceTerm(otherSide,iFace);

          // reconstruct the current face gradients
          m_faceTermComputer->
              reconstructFluxPntsGradients(otherSide,m_pertGrads[otherSide],m_grads[otherSide].size());

          // compute the perturbed face term
          m_faceTermComputer->computeDiffFaceTerm(m_pertResUpdates);

          // compute the finite difference derivative of the face term
          m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

          // multiply residual update derivatives with residual factor and +1.0 or -1.0 depending on otherSide
          m_derivResUpdates *= pertSideFac;

          // add the derivative of the residual updates to the accumulator
          CFuint resUpdIdx = 0;
          for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace, resUpdIdx += m_nbrEqs)
          {
            CFuint cvIdx = (*m_cvCVConnSVFace)[m_orient][iSubFace][side];
            accSC.addValues(cvIdx,iCVPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
          }

          // restore physical variable in state
          m_numJacob->restore(pertState[iEqPert]);

          // restore physical variable in the flux points (contributing face term computer)
          contrFaceTermComputer->restorePhysVar(contrCellSide,iEqPert);
        }
      }
// accSC.printToScreen();

      if (getMethodData().doComputeJacobian())
      {
        // add the values to the jacobian matrix
        m_lss->getMatrix()->addValues(accSC);
       }

      // reset to zero the entries in the block accumulator
      accSC.reset();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::setup()
{
  CFAUTOTRACE;

  DiffFaceTermRHSSpectralFV::setup();

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

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get number of equations and dimensionality
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  cf_assert(svLocalData.size() > 0);

  // get the number of CVs in a cell
  const CFuint nbrCVs = svLocalData[0]->getNbrOfCVs();

  // maximum number of CVs at SV face
  CFuint maxNbrCVsAtSVFace = 0;
  SafePtr< vector< vector< CFuint > > > svFaceCVConn = svLocalData[0]->getExtSVFaceCVConn();
  for (CFuint iFace = 0; iFace < svFaceCVConn->size(); ++iFace)
  {
    const CFuint nbrCVsAtSVFace = (*svFaceCVConn)[iFace].size();
    maxNbrCVsAtSVFace = maxNbrCVsAtSVFace > nbrCVsAtSVFace ? maxNbrCVsAtSVFace : nbrCVsAtSVFace;
  }

  // get the number of faces in a cell
  const CFuint nbrFaces = svLocalData[0]->getNbrSVFaces();
  const CFuint nbrFacesM1 = nbrFaces - 1;

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(2*nbrCVs,2*nbrCVs,m_nbrEqs));

  // create single cell blockaccumulator
  m_accSC.reset(m_lss->createBlockAccumulator(nbrCVs,nbrCVs,m_nbrEqs));

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
  const CFuint nbrFaceFluxes = maxNbrCVsAtSVFace*m_nbrEqs;
  m_pertResUpdates .resize(nbrFaceFluxes);
  m_derivResUpdates.resize(nbrFaceFluxes);
  const CFuint nbrCellResiduals = nbrCVs*m_nbrEqs;
  m_unpertCellDiffRes       .resize(2);
  m_unpertCellDiffRes[LEFT ].resize(nbrCellResiduals);
  m_unpertCellDiffRes[RIGHT].resize(nbrCellResiduals);
  m_pertCellDiffRes         .resize(nbrCellResiduals);
  m_derivCellDiffRes        .resize(nbrCellResiduals);

  // allocate memory for perturbed gradiens
  m_pertGrads.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_pertGrads[iSide].resize(nbrCVs);
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
    {
      m_pertGrads[iSide][iCV] = new vector<RealVector>(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        (*m_pertGrads[iSide][iCV])[iEq].resize(dim);
      }
    }
  }

  // resize m_cellGradsMinusFaceTerm
  m_cellGradsMinusFaceTerm.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_cellGradsMinusFaceTerm[iSide].resize(nbrCVs);
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
    {
      m_cellGradsMinusFaceTerm[iSide][iCV].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_cellGradsMinusFaceTerm[iSide][iCV][iEq].resize(dim);
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
      m_cellGradsMinusOtherFaceTerm[iSide][iFace].resize(nbrCVs);
      for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
      {
        m_cellGradsMinusOtherFaceTerm[iSide][iFace][iCV].resize(m_nbrEqs);
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_cellGradsMinusOtherFaceTerm[iSide][iFace][iCV][iEq].resize(dim);
        }
      }
    }
  }

  // resize gradient updates
  m_gradUpdates.resize(nbrCVs);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_gradUpdates[iCV].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradUpdates[iCV][iEq].resize(dim);
    }
  }

  // resize neighbouring SVs CV volumes
  m_cvInvVols.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_cvInvVols[iSide].resize(nbrCVs);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSJacobSpectralFV::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSide = 0; iSide < m_pertGrads.size(); ++iSide)
  {
    for (CFuint iCV = 0; iCV < m_pertGrads[iSide].size(); ++iCV)
    {
      deletePtr(m_pertGrads[iSide][iCV]);
    }
    m_pertGrads[iSide].resize(0);
  }
  m_pertGrads.resize(0);

  DiffFaceTermRHSSpectralFV::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
