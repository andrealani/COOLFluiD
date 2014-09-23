#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDLinEuler/LinEulerMeanFlowSourceTermRHSSpectralFD.hh"
#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::LinearizedEuler;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

    MethodCommandProvider<LinEulerMeanFlowSourceTermRHSSpectralFD, SpectralFDMethodData, SpectralFDLinEulerModule> LinEulerMeanFlowSourceTermRHSSpectralFDProvider("LinEulerMeanFlow");

//////////////////////////////////////////////////////////////////////////////

LinEulerMeanFlowSourceTermRHSSpectralFD::LinEulerMeanFlowSourceTermRHSSpectralFD(const std::string& name) :
  StdSourceTerm(name),
  m_linEulerVarSet(CFNULL),
  m_cellBuilder(CFNULL),
  m_volTermComputer(CFNULL),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_bcStateComputers(),
  m_iElemType(),
  m_cellExtraVarGrads(),
  m_faceNghbrStates(),
  m_isFaceOnBoundary(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_gradUpdates(),
  m_dim(),
  m_srcTerm()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

LinEulerMeanFlowSourceTermRHSSpectralFD::~LinEulerMeanFlowSourceTermRHSSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::configure ( Config::ConfigArgs& args )
{
  StdSourceTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::getSourceTermData()
{
  // call parent class function
  StdSourceTerm::getSourceTermData();

  // set data in volume and face term computers
  setVolumeAndFaceTermComputersData();

  // resize gradient updates
  resizeGradUpdates();

  // get the geodata of the geometric entity builder and set the TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = trs;
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::addSourceTerm()
{
  // get the geodata of the geometric entity builder and set the element index
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.idx = m_cell->getID();

  // rebuild the cell with Cells2Faces builder
  m_cell = m_cellBuilder->buildGE();

  // COMPUTE DATA IN VOLUME AND FACE TERM COMPUTERS
  computeCellVolumeAndFaceTermData();

  // COMPUTE EXTRA VARIABLE GRADIENTS
  computeExtraVarsGradients();

  // GET EXTRA VARIABLES AT THE SOLUTION POINTS
//  const vector< RealVector* >& cellExtraVars = *m_volTermComputer->getCellExtraVars();

  // COMPUTE SOURCE TERMS
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // loop over solution points in this cell to add the source term
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrSol = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
/*============================================
    // get the source term from the linearized Euler variable set
    m_linEulerVarSet->computeMeanFlowSourceTerm(*(*m_cellStates)[iSol],
                                                *cellExtraVars[iSol],
                                                *m_cellExtraVarGrads[iSol],
                                                m_srcTerm);
//============================================*/
    // loop over physical variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
    {
      rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
    }
  }

  //release the GeometricEntity
  m_cellBuilder->releaseGE();
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::setVolumeAndFaceTermComputersData()
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
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::resizeGradUpdates()
{
  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

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

void LinEulerMeanFlowSourceTermRHSSpectralFD::computeCellVolumeAndFaceTermData()
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
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::computeExtraVarsGradients()
{
  // compute volume contribution to current cell extra variable gradients
  m_volTermComputer->computeGradientExtraVarsVolumeTerm(m_gradUpdates[0]);
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellExtraVarGrads[iSol])[iEq] = m_gradUpdates[0][iSol][iEq];
    }
  }

  // compute face contributions to the gradients
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // compute the boundary face contribution to the gradients of the extra variables
      m_bndFaceTermComputers[iFace]->computeGradientExtraVarsFaceTerm(m_gradUpdates[0]);

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellExtraVarGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
        }
      }
    }
    else
    {
      // compute the face contribution to the gradients of the extra variables
      m_faceTermComputers[iFace]->computeGradientExtraVarsFaceTerm(m_gradUpdates);

      // add the contribution to the gradients
      const CFuint side = (*m_currCellSide)[iFace];
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellExtraVarGrads[iSol])[iEq] += m_gradUpdates[side][iSol][iEq];
        }
      }
    }
  }

  // divide by solution point Jacobian determinant
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solPntJacobDets[iSol];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellExtraVarGrads[iSol])[iEq] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  StdSourceTerm::setup();

  // get linearized Euler variable set
  m_linEulerVarSet = getMethodData().getUpdateVar().d_castTo<LinEulerVarSet>();

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

  // get the dimensionality
  m_dim    = PhysicalModelStack::getActive()->getDim();

  // get cell builder
  m_cellBuilder = getMethodData().getCellBuilder();
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
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
    m_faceTermComputers   [iFace] =
        getMethodData().getAdditionalFaceTermComputer   (faceCompIdx);
    m_bndFaceTermComputers[iFace] =
        getMethodData().getAdditionalBndFaceTermComputer(faceCompIdx);
  }

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(maxNbrFaces);
  for (CFuint iFace = 0; iFace < maxNbrFaces; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

    // resize m_gradUpdates
  m_gradUpdates.resize(2);

  // allocate memory for m_cellExtraVarGrads
  m_cellExtraVarGrads.resize(maxNbrSolPnts);
  for (CFuint iSol = 0; iSol < maxNbrSolPnts; ++iSol)
  {
    m_cellExtraVarGrads[iSol] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellExtraVarGrads[iSol])[iEq].resize(m_dim);
    }
  }

  // resize m_srcTerm
  m_srcTerm.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void LinEulerMeanFlowSourceTermRHSSpectralFD::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSol = 0; iSol < m_cellExtraVarGrads.size(); ++iSol)
  {
    deletePtr(m_cellExtraVarGrads[iSol]);
  }
  m_cellExtraVarGrads.resize(0);

  // call unsetup of parent class
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
