#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/VolTermRHSSpectralFD.hh"

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

MethodCommandProvider<VolTermRHSSpectralFD, SpectralFDMethodData, SpectralFDModule> VolTermRHSSpectralFDProvider("VolTermRHS");

//////////////////////////////////////////////////////////////////////////////

VolTermRHSSpectralFD::VolTermRHSSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhs("rhs"),
  m_cellBuilder(CFNULL),
  m_volTermComputer(CFNULL),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_bcStateComputers(),
  m_iElemType(),
  m_cell(CFNULL),
  m_cellStates(),
  m_cellGrads(),
  m_faceNghbrStates(),
  m_isFaceOnBoundary(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_resUpdates(),
  m_diffResUpdates(),
  m_gradUpdates(),
  m_bndFaceGradUpdates(),
  m_solJacobDet(),
  m_solPntsLocalCoords(CFNULL),
  m_nbrEqs(),
  m_dim()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

VolTermRHSSpectralFD::~VolTermRHSSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // SET THE VOLUME AND FACE TERM COMPUTERS DATA FOR THIS ELEMENT TYPE
    setVolumeAndFaceTermComputersData();

    // RESIZE THE VARIABLES M_RESUPDATES AND M_GRADUPDATES
    resizeResAndGradUpdates();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // COMPUTE DATA IN VOLUME AND FACE TERM COMPUTERS
        computeCellVolumeAndFaceTermData();

        // COMPUTE CONVECTIVE VOLUME TERMS
        m_volTermComputer->computeCellConvVolumeTerm(m_resUpdates);

        // COMPUTE CELL GRADIENTS
        computeAndReconstructGradients();

        // COMPUTE DIFFUSIVE VOLUME TERMS
        m_volTermComputer->computeCellDiffVolumeTerm(m_diffResUpdates);

        // ADD TOTAL UPDATE TO RESIDUAL
        m_resUpdates += m_diffResUpdates;
        addUpdatesToResidual();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  CFTRACEEND;
 }

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::setVolumeAndFaceTermComputersData()
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

void VolTermRHSSpectralFD::resizeResAndGradUpdates()
{
  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // number of entries in residual updates
  const CFuint nbrResUpdates = m_nbrEqs*nbrSolPnts;

  // resize m_resUpdates and m_diffResUpdates
  m_resUpdates    .resize(nbrResUpdates);
  m_diffResUpdates.resize(nbrResUpdates);

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

  // resize m_bndFaceGradUpdates
  const CFuint maxNbrBndFaces = m_bndFaceGradUpdates.size();
  for (CFuint iBnd = 0; iBnd < maxNbrBndFaces; ++iBnd)
  {
    m_bndFaceGradUpdates[iBnd].resize(nbrSolPnts);
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_bndFaceGradUpdates[iBnd][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_bndFaceGradUpdates[iBnd][iSol][iEq].resize(m_dim);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::computeCellVolumeAndFaceTermData()
{
  // CELL DATA IN VOLUME TERM COMPUTER
  // set the cell in the volume term computer
  m_volTermComputer->setCurrentCell(m_cell);

  // compute the cell data in the volume term computer
  m_volTermComputer->computeCellData();

  // reconstruct the flux point states in the volume term computer
  m_volTermComputer->reconstructStates(*m_cellStates);

  // FACE AND CELL DATA IN FACE TERM COMPUTERS
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  const vector<GeometricEntity*>* faces = m_cell->getNeighborGeos();
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
    }
  }

  // compute solution points Jacobian determinants
  m_solJacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::addUpdatesToResidual()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // update rhs
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrRes = m_resUpdates.size();
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes, ++resID)
  {
      rhs[resID] += resFactor*m_resUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::computeAndReconstructGradients()
{
  // compute volume contribution to current cell gradients
  m_volTermComputer->computeGradientVolumeTerm(m_gradUpdates[0]);
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[iSol])[iEq] = m_gradUpdates[0][iSol][iEq];
    }
  }

  // compute face contributions to the gradients
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // compute the boundary face contribution to the gradients
      m_bndFaceTermComputers[iFace]->computeGradientFaceTerm(m_bndFaceGradUpdates[iFace]);

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[iSol])[iEq] += m_bndFaceGradUpdates[iFace][iSol][iEq];
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
          (*m_cellGrads[iSol])[iEq] += m_gradUpdates[side][iSol][iEq];
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
      (*m_cellGrads[iSol])[iEq] *= invJacobDet;
    }
  }

  // reconstruct the gradients in the flux points
  m_volTermComputer->reconstructGradients(m_cellGrads);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::setup()
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

  // get the number of equations in the physical modes and the dimensionality
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim();

  // get cell builder
  m_cellBuilder = getMethodData().getCellBuilder();;
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_currCellSide     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrients      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdx        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  // get the volume term computer
  m_volTermComputer = getMethodData().getVolTermComputer().d_castTo<CompactVolTermComputer>();

  // get face term computers and boundary face term computers
  m_faceTermComputers   .resize(maxNbrFaces);
  m_bndFaceTermComputers.resize(maxNbrFaces);
  m_faceTermComputers   [0] = getMethodData().getFaceTermComputer   ().d_castTo<CompactFaceTermComputer>();
  m_bndFaceTermComputers[0] = getMethodData().getBndFaceTermComputer().d_castTo<CompactBndFaceTermComputer>();
  CFuint faceCompIdx = 0;
  for (CFuint iFace = 1; iFace < maxNbrFaces; ++iFace, ++faceCompIdx)
  {
    m_faceTermComputers   [iFace] =
        getMethodData().getAdditionalFaceTermComputer   (faceCompIdx).d_castTo<CompactFaceTermComputer>();
    m_bndFaceTermComputers[iFace] =
        getMethodData().getAdditionalBndFaceTermComputer(faceCompIdx).d_castTo<CompactBndFaceTermComputer>();
  }

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

    // resize m_gradUpdates and m_bndFaceGradUpdates
  m_gradUpdates.resize(2);
  m_bndFaceGradUpdates.resize(maxNbrFaces);

    // resize m_solJacobDet
  m_solJacobDet.resize(maxNbrSolPnts);

  // allocate memory for m_currCellGrads
  m_cellGrads.resize(maxNbrSolPnts);
  for (CFuint iSol = 0; iSol < maxNbrSolPnts; ++iSol)
  {
    m_cellGrads[iSol] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[iSol])[iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSSpectralFD::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSol = 0; iSol < m_cellGrads.size(); ++iSol)
  {
    deletePtr(m_cellGrads[iSol]);
  }
  m_cellGrads.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
VolTermRHSSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
