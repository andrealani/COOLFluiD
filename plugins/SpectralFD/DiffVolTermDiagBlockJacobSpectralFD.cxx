#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/DiffVolTermDiagBlockJacobSpectralFD.hh"
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

MethodCommandProvider<DiffVolTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule>
    DiffVolTermDiagBlockJacobSpectralFDProvider("DiffVolTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

DiffVolTermDiagBlockJacobSpectralFD::DiffVolTermDiagBlockJacobSpectralFD(const std::string& name) :
  DiffVolTermRHSJacobSpectralFD(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_currDiagMatrix(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffVolTermDiagBlockJacobSpectralFD::~DiffVolTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermDiagBlockJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

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
        // set current diagonal block Jacobian matrix
        m_currDiagMatrix = &diagBlockJacobMatr[elemIdx];

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

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermDiagBlockJacobSpectralFD::setVolumeTermData()
{
  // call the setVolumeTermData() of the parent class of the parent class,
  // such that the accumulator for the lss matrix is not set.
  DiffVolTermRHSSpectralFD::setVolumeTermData();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get solution point local coordinates
  m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();

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

void DiffVolTermDiagBlockJacobSpectralFD::computeJacobDiffVolTerm()
{
  // get number of solution points in this cell
  const CFuint nbrSolPnts = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // loop over the states/solpnts in this cell to perturb the states
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

      // add the derivative of the residual updates to diagonal block matrix
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points
      m_volTermComputer->restorePhysVar(iEqPert);

      // restore physical variable in the face flux points
      restoreFacePhysVars(iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes == m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class of parent class, in order not to call the lss
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

void DiffVolTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSol = 0; iSol < m_pertGrads.size(); ++iSol)
  {
    deletePtr(m_pertGrads[iSol]);
  }
  m_pertGrads.resize(0);


  // call unsetup of parent class of parent class, in order not to call the lss
  DiffVolTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffVolTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffVolTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr     );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
