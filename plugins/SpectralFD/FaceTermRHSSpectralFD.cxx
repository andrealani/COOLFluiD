#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/FaceTermRHSSpectralFD.hh"
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

MethodCommandProvider<FaceTermRHSSpectralFD, SpectralFDMethodData, SpectralFDModule> FaceTermRHSSpectralFDProvider("FaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

FaceTermRHSSpectralFD::FaceTermRHSSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  m_faceTermComputer(CFNULL),
  m_faceBuilder(CFNULL),
  m_face(),
  m_cells(),
  m_orient(),
  m_states(),
  m_resUpdates(),
  m_diffResUpdates(),
  m_updateCoefContr(),
  m_diffUpdateCoefContr(),
  m_nbrEqs()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

FaceTermRHSSpectralFD::~FaceTermRHSSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::execute()
{
  CFTRACEBEGIN;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

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

  // SET THE DATA NEEDED TO COMPUTE THE FACE TERMS
  setFaceTermComputerData();

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // SET THE ORIENTATION OF THE FACES
    m_faceTermComputer->setFaceOrientation(m_orient);

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoData.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      // if one of the neighbouring cells is parallel updatable,
      // compute the face terms
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // COMPUTE DATA IN FACE TERM COMPUTER
        computeFaceTermData();

        // COMPUTE CONVECTIVE FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
        m_faceTermComputer->
            computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_updateCoefContr);

        // COMPUTE DIFFUSIVE FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
        m_faceTermComputer->
            computeDiffFaceTermAndUpdateCoefContributions(m_diffResUpdates,m_diffUpdateCoefContr);

        // ADD TOTAL UPDATE TO RESIDUAL AND UPDATE COEFFICIENTS
        m_resUpdates[LEFT ] += m_diffResUpdates[LEFT ];
        m_resUpdates[RIGHT] += m_diffResUpdates[RIGHT];
        addUpdatesToResidual();
        m_updateCoefContr[LEFT ] += m_diffUpdateCoefContr[LEFT ];
        m_updateCoefContr[RIGHT] += m_diffUpdateCoefContr[RIGHT];
        addUpdateCoeffContributions();
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
// CF_DEBUG_EXIT;
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::setFaceTermComputerData()
{
  // set the face term data in the face term computer
  m_faceTermComputer->setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::computeFaceTermData()
{
  // set the face in the boundary face term computer
  m_faceTermComputer->setCurrentFace(m_face);

  // compute the face data in the face term computer
  m_faceTermComputer->computeFaceData();

  // compute the neighbouring cell data in the face term computer
  m_faceTermComputer->computeNeighbourCellData();

  // reconstruct the states in the face flux points
  m_faceTermComputer->reconstructFluxPntsStates(m_states);

  // compute the solution polynomial gradients in the face flux points
  m_faceTermComputer->reconstructFluxPntsSolPolyGrads(m_states);
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::addUpdatesToResidual()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // add updates to residuals
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    CFuint resID = m_nbrEqs*( (*m_states[iSide])[0]->getLocalID() );

    const CFuint nbrRes = m_resUpdates[iSide].size();
    /// @warning this size is for the highest polynomial degree, this should probably be changed for p-multigrid
    for (CFuint iRes = 0; iRes < nbrRes; ++iRes, ++resID)
    {
      rhs[resID] += resFactor*m_resUpdates[iSide][iRes];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::addUpdateCoeffContributions()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();
      updateCoeff[solID] += m_updateCoefContr[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::setup()
{
  CFAUTOTRACE;

  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_faceTermComputer = getMethodData().getFaceTermComputer().d_castTo<CompactFaceTermComputer>();

  // number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // maximum number of solution points in a cell
  const CFuint maxNbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  m_resUpdates.resize(2);
  m_resUpdates[LEFT ].resize(maxNbrSolPnts*m_nbrEqs);
  m_resUpdates[RIGHT].resize(maxNbrSolPnts*m_nbrEqs);
  m_diffResUpdates.resize(2);
  m_diffResUpdates[LEFT ].resize(maxNbrSolPnts*m_nbrEqs);
  m_diffResUpdates[RIGHT].resize(maxNbrSolPnts*m_nbrEqs);
  m_states.resize(2);
  m_updateCoefContr.resize(2);
  m_diffUpdateCoefContr.resize(2);
  m_cells.resize(2);
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermRHSSpectralFD::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FaceTermRHSSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
