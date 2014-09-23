#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/DiffVolTermRHSJacobSpectralFD.hh"
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

MethodCommandProvider<DiffVolTermRHSJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> DiffVolTermRHSJacobSpectralFDProvider("DiffVolTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

DiffVolTermRHSJacobSpectralFD::DiffVolTermRHSJacobSpectralFD(const std::string& name) :
  DiffVolTermRHSSpectralFD(name),
  m_cellBuilder(CFNULL),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(),
  m_faces(CFNULL),
  m_faceNghbrStates(),
  m_pertResUpdates(),
  m_derivResUpdates(),
  m_gradUpdates(),
  m_pertGrads(),
  m_solJacobDet(),
  m_solPntsLocalCoords(CFNULL),
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

DiffVolTermRHSJacobSpectralFD::~DiffVolTermRHSJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  DiffVolTermRHSSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFD::execute()
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
      // build the cell
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the neighbouring faces
      m_faces = m_cell->getNeighborGeos();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // get the gradients in this cell (in variable m_cellGrads)
        setGradients();

        // set the current cell and compute the cell data in the volume term computer
        m_volTermComputer->setCurrentCell(m_cell);
        m_volTermComputer->computeCellData();

        // reconstruct the solution in the flux points
        m_volTermComputer->reconstructStates(*m_cellStates);

        // reconstruct the gradients in the flux points
        m_volTermComputer->reconstructGradients(m_cellGrads);

        // compute the volume term
        m_volTermComputer->computeCellDiffVolumeTerm(m_resUpdates);

        // update rhs
        updateRHS();

        // get the neigbouring states of the faces
        setFaceNeighbourStates();

        // compute the data of the cell
        setCellData();

        // compute the contribution to the jacobian
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

void DiffVolTermRHSJacobSpectralFD::setVolumeTermData()
{
  DiffVolTermRHSSpectralFD::setVolumeTermData();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get solution point local coordinates
  m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(nbrSolPnts,nbrSolPnts,m_nbrEqs));

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

void DiffVolTermRHSJacobSpectralFD::resizeResAndGradUpdates()
{
  // call function in parent class
  DiffVolTermRHSSpectralFD::resizeResUpdatesAndCellGrads();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // number of entries in residual updates
  const CFuint nbrResUpdates = m_nbrEqs*nbrSolPnts;

  // resize m_pertResUpdates and m_derivResUpdate
  m_pertResUpdates .resize(nbrResUpdates);
  m_derivResUpdates.resize(nbrResUpdates);

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

void DiffVolTermRHSJacobSpectralFD::setFaceNeighbourStates()
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

void DiffVolTermRHSJacobSpectralFD::setCellData()
{
  // compute solution points Jacobian determinants
  m_solJacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // set orientation or boundary condition state computers of the neighbouring faces
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
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

      // compute the neighbouring cell data in the face term computer
//       m_faceTermComputers[iFace]->computeNeighbourCellData();

      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->reconstructFluxPntsStates(m_faceNghbrStates[iFace]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFD::backupAndReconstructFacePhysVars(const CFuint iVar)
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

void DiffVolTermRHSJacobSpectralFD::restoreFacePhysVars(const CFuint iVar)
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

void DiffVolTermRHSJacobSpectralFD::computePerturbedGradients()
{
  // compute the internal contribution to the gradients
  m_volTermComputer->computeGradientVolumeTerm(m_gradUpdates[0]);
  const CFuint nbrSolPnts = m_cellGrads.size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[iSol])[iEq] = m_gradUpdates[0][iSol][iEq];
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
          (*m_pertGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
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
          (*m_pertGrads[iSol])[iEq] += m_gradUpdates[side][iSol][iEq];
        }
      }
    }
  }

  // divide by solution point Jacobian determinant
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDet[iSol];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[iSol])[iEq] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFD::computeJacobDiffVolTerm()
{
  // get number of solution points in this cell
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

  // loop over the states in this cell to perturb the states
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
      m_volTermComputer->backupAndReconstructPhysVar(iEqPert,*m_cellStates);

      // backup and reconstruct physical variable in the face flux points
      backupAndReconstructFacePhysVars(iEqPert);

      // compute the perturbed cell gradients after the perturbation
      computePerturbedGradients();

      // reconstruct the gradients in the flux points
      m_volTermComputer->reconstructGradients(m_pertGrads);

      // compute the perturbed diffusive volume term
      m_volTermComputer->computeCellDiffVolumeTerm(m_pertResUpdates);

      // compute the finite difference derivative of the volume term
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

void DiffVolTermRHSJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  DiffVolTermRHSSpectralFD::setup();

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

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim();

  // get face term computers and boundary face term computers
  m_faceTermComputers   .resize(maxNbrFaces);
  m_bndFaceTermComputers.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceTermComputers   [iFace] = getMethodData().getAdditionalFaceTermComputer   (iFace);
    m_bndFaceTermComputers[iFace] = getMethodData().getAdditionalBndFaceTermComputer(iFace);
  }

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

  // resize m_gradUpdates
  m_gradUpdates.resize(2);

  // resize m_solJacobDet
  m_solJacobDet.resize(maxNbrSolPnts);

  // allocate memory for m_pertGrads
  m_pertGrads.resize(maxNbrSolPnts);
  for (CFuint iSol = 0; iSol < maxNbrSolPnts; ++iSol)
  {
    m_pertGrads[iSol] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_pertGrads[iSol])[iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSol = 0; iSol < m_pertGrads.size(); ++iSol)
  {
    deletePtr(m_pertGrads[iSol]);
  }
  m_pertGrads.resize(0);

  DiffVolTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
