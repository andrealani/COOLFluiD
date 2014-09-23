#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFV/DiffVolTermRHSJacobSpectralFV.hh"
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

MethodCommandProvider<DiffVolTermRHSJacobSpectralFV, SpectralFVMethodData, SpectralFVModule> DiffVolTermRHSJacobSpectralFVProvider("DiffVolTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

DiffVolTermRHSJacobSpectralFV::DiffVolTermRHSJacobSpectralFV(const std::string& name) :
  DiffVolTermRHSSpectralFV(name),
  m_cellBuilder(CFNULL),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(),
  m_faces(),
  m_faceNghbrStates(),
  m_pertResUpdates(),
  m_derivResUpdates(),
  m_gradUpdates(),
  m_pertGrads(),
  m_cvInvVols(),
  m_invVolFracCVs(),
  m_cvCVConnSVFace(),
  m_svFaceCVConn(),
  m_isFaceOnBoundary(CFNULL),
  m_nghbrCellSide(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_bcStateComputers(CFNULL),
  m_dim()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

DiffVolTermRHSJacobSpectralFV::~DiffVolTermRHSJacobSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::configure ( Config::ConfigArgs& args )
{
  DiffVolTermRHSSpectralFV::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the cell builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // set the volume term data for this element type
    setVolumeTermData();

    // resize the variables m_resUpdates and m_cellGrads
    resizeResAndGradUpdates();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the neighbouring faces
      m_faces = m_cell->getNeighborGeos();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // get the neigbouring states of the faces
        setFaceNeighbourStates();

        // get the gradients in this cell (in variable m_cellGrads)
        setGradients();

        // set the current cell and compute the cell data in the volume term computer
        m_volTermComputer->setCurrentCell(m_cell);
        m_volTermComputer->computeCellData();

        // compute the data of the cell
        setCellData();

        // reconstruct the solution in the flux points
        m_volTermComputer->reconstructStates(*m_cellStates);

        // reconstruct the gradients in the flux points
        m_volTermComputer->reconstructGradients(m_cellGrads,m_cellGrads.size());

        // compute the volume term
        m_volTermComputer->computeCellDiffVolumeTerm(m_resUpdates);

        // update rhs
        updateRHS();

        // compute the diffusive volume term contribution to the jacobian
        computeJacobDiffVolTerm();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
// CF_DEBUG_EXIT;
  CFTRACEEND;
 }

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::setVolumeTermData()
{
  DiffVolTermRHSSpectralFV::setVolumeTermData();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get the number of CVs in current element type
  const CFuint nbrCVs = svLocalData[m_iElemType]->getNbrOfCVs();

  // get the volume fractions of the CVs
  m_invVolFracCVs = svLocalData[m_iElemType]->getInvVolFracCV();

  // get the SV faces-CV connectivity
  m_cvCVConnSVFace = svLocalData[m_iElemType]->getExtSVFaceCVConnPerOrient();

  // get the boundary SV face-CV connectivity
  m_svFaceCVConn = svLocalData[m_iElemType]->getExtSVFaceCVConn();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(nbrCVs,nbrCVs,m_nbrEqs));

  // set the data in the face term and boundary face term computers
  for (CFuint iFace = 0; iFace < m_faceTermComputers.size(); ++iFace)
  {
    m_faceTermComputers[iFace]->setFaceTermData();
  }
  for (CFuint iFace = 0; iFace < m_bndFaceTermComputers.size(); ++iFace)
  {
    m_bndFaceTermComputers[iFace]->setFaceTermData();
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::resizeResAndGradUpdates()
{
  // call function in parent class
  DiffVolTermRHSSpectralFV::resizeResUpdatesAndCellGrads();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // number of CVs
  const CFuint nbrCVs = svLocalData[m_iElemType]->getNbrOfCVs();

  // number of entries in residual updates
  const CFuint nbrResUpdates = m_nbrEqs*nbrCVs;

  // resize m_pertResUpdates and m_derivResUpdate
  m_pertResUpdates .resize(nbrResUpdates);
  m_derivResUpdates.resize(nbrResUpdates);

  // resize m_gradUpdates
  m_gradUpdates.resize(nbrCVs);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_gradUpdates[iCV].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradUpdates[iCV][iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::setFaceNeighbourStates()
{
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    const CFuint nbrFaceNeighbours = (*m_faces)[iFace]->nbNeighborGeos();
    for (CFuint iSide = 0; iSide < nbrFaceNeighbours; ++iSide)
    {
      GeometricEntity* cell = (*m_faces)[iFace]->getNeighborGeo(iSide);
      m_faceNghbrStates[iFace][iSide] = cell->getStates();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::setCellData()
{
  // cell volume
  const CFreal cellInvVol = 1.0/m_cell->computeVolume();

  // compute CV volumes
  const CFuint nbrCVs = m_invVolFracCVs->size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_cvInvVols[iCV] = cellInvVol*(*m_invVolFracCVs)[iCV];
  }

  // set orientation or boundary condition state computers of the neighbouring faces
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // set BC state computer
      m_bndFaceTermComputers[iFace]->setBcStateComputer((*m_bcStateComputers)[(*m_faceBCIdx)[iFace]]);

      // set the orientation of the face (set during setup)

      // set the face in the boundary face term computer
      m_bndFaceTermComputers[iFace]->setCurrentFace((*m_faces)[iFace]);

      // compute the face data in the boundary face term computer
      m_bndFaceTermComputers[iFace]->computeFaceData();

//       // compute the neighbouring cell data in the boundary face term computer
//       m_bndFaceTermComputers[iFace]->computeNeighbourCellData();

      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->reconstructFluxPntsStates(*m_cellStates);
    }
    else
    {
      // set the orientation of the face
      m_faceTermComputers[iFace]->setFaceOrientation((*m_faceOrients)[iFace]);

      // set the face in the face term computer
      m_faceTermComputers[iFace]->setCurrentFace((*m_faces)[iFace]);

      // compute the face data in the face term computer
      m_faceTermComputers[iFace]->computeFaceData();

//       // compute the neighbouring cell data in the face term computer
//       m_faceTermComputers[iFace]->computeNeighbourCellData();

      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->reconstructFluxPntsStates(*m_faceNghbrStates[iFace][LEFT ],
                                                            *m_faceNghbrStates[iFace][RIGHT]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::backupAndReconstructFacePhysVars(const CFuint iVar)
{
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->backupAndReconstructPhysVar(iVar,*m_cellStates);
    }
    else
    {
      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->backupAndReconstructPhysVar((*m_currCellSide)[iFace],iVar,
                                                              *m_cellStates);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::restoreFacePhysVars(const CFuint iVar)
{
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->restorePhysVar(iVar);
    }
    else
    {
      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->restorePhysVar((*m_currCellSide)[iFace],iVar);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::computePerturbedGradients()
{
  // compute the internal contribution to the gradients
  m_volTermComputer->computeGradientVolumeTerm(m_gradUpdates);
  const CFuint nbrCVs = m_cellGrads.size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[iCV])[iEq] = m_gradUpdates[iCV][iEq];
    }
  }

  // compute face contributions to the gradients
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // compute the boundary face contribution to the gradients
      m_bndFaceTermComputers[iFace]->computeGradientFaceTerm(m_gradUpdates);

      // add the contribution to the gradients
      const CFuint nbrSubFaces = (*m_svFaceCVConn)[iFace].size();
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
      {
        // get internal face neighbour CV index
        const CFuint cvIdx  = (*m_svFaceCVConn)[iFace][iSubFace];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_pertGrads[cvIdx])[iEq] += m_gradUpdates[iSubFace][iEq];
        }
      }
    }
    else
    {
      // compute the boundary face contribution to the gradients
      m_faceTermComputers[iFace]->computeGradientFaceTerm(m_gradUpdates);

      // add the contribution to the gradients
      const CFreal sideFactor = pow(-1.0,static_cast<CFreal>((*m_currCellSide)[iFace]));
      const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[(*m_faceOrients)[iFace]].size();
      for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
      {
        // get internal face neighbour CV index
        const CFuint cvIdx  = (*m_cvCVConnSVFace)[(*m_faceOrients)[iFace]][iSubFace][(*m_currCellSide)[iFace]];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_pertGrads[cvIdx])[iEq] += sideFactor*m_gradUpdates[iSubFace][iEq];
        }
      }
    }
  }

  // divide by CV volumes
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[iCV])[iEq] *= m_cvInvVols[iCV];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::computeJacobDiffVolTerm()
{
  // get number of CVs in this cell
  const CFuint nbrCVs = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    acc.setRowColIndex(iCV,(*m_cellStates)[iCV]->getLocalID());
  }

  // loop over the states/CVs in this cell to perturb the states
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
      m_volTermComputer->backupAndReconstructPhysVar(iEqPert,*m_cellStates);

      // backup and reconstruct physical variable in the face flux points
      backupAndReconstructFacePhysVars(iEqPert);

      // compute the perturbed cell gradients after the perturbation
      computePerturbedGradients();

      // reconstruct the gradients in the flux points
      m_volTermComputer->reconstructGradients(m_pertGrads,m_cellGrads.size());

      // compute the perturbed diffusive volume term
      m_volTermComputer->computeCellDiffVolumeTerm(m_pertResUpdates);

      // compute the finite difference derivative of the volume term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iCV = 0; iCV < nbrCVs; ++iCV, resUpdIdx += m_nbrEqs)
      {
        acc.addValues(iCV,iCVPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points
      m_volTermComputer->restorePhysVar(iEqPert);

      // restore physical variable in the face flux points
      restoreFacePhysVars(iEqPert);
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
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::setup()
{
  CFAUTOTRACE;

  DiffVolTermRHSSpectralFV::setup();

  // get CellToFaceGeBuilder
  m_cellBuilder      = getMethodData().getCellBuilder();
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSide    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSide     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrients      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdx        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // maximum number of CVs and faces in a cell
  CFuint maxNbrCVs = 0;
  CFuint maxNbrFaces = 0;
  for (CFuint iElemType = 0; iElemType < svLocalData.size(); ++iElemType)
  {
    const CFuint nbrCVs = svLocalData[iElemType]->getNbrOfCVs();
    maxNbrCVs = maxNbrCVs > nbrCVs ? maxNbrCVs : nbrCVs;

    const CFuint nbrFaces = svLocalData[iElemType]->getNbrSVFaces();
    maxNbrFaces = maxNbrFaces > nbrFaces ? maxNbrFaces : nbrFaces;
  }

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim();

  // get face term computers and boundary face term computers
  m_faceTermComputers   .resize(maxNbrFaces);
  m_bndFaceTermComputers.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceTermComputers   [iFace] = getMethodData().getAdditionalFaceTermComputer   (iFace);
    m_bndFaceTermComputers[iFace] = getMethodData().getAdditionalBndFaceTermComputer(iFace);
    m_bndFaceTermComputers[iFace]->setFaceOrientation(iFace);
  }

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

  // resize m_cvInvVols
  m_cvInvVols.resize(maxNbrCVs);

  // allocate memory for m_pertGrads
  m_pertGrads.resize(maxNbrCVs);
  for (CFuint iCV = 0; iCV < maxNbrCVs; ++iCV)
  {
    m_pertGrads[iCV] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[iCV])[iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFV::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iCV = 0; iCV < m_pertGrads.size(); ++iCV)
  {
    deletePtr(m_pertGrads[iCV]);
  }
  m_pertGrads.resize(0);

  DiffVolTermRHSSpectralFV::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
