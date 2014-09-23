#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/ConvFaceTermRHSSpectralFD.hh"
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

MethodCommandProvider<ConvFaceTermRHSSpectralFD, SpectralFDMethodData, SpectralFDModule> ConvFaceTermRHSSpectralFDProvider("ConvFaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermRHSSpectralFD::ConvFaceTermRHSSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  m_faceBuilder(CFNULL),
  m_faceTermComputer(CFNULL),
  m_face(),
  m_cells(),
  m_orient(),
  m_states(),
  m_resUpdates(),
  m_waveSpeedUpd(),
  m_gradUpdates(),
  m_nbrEqs()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermRHSSpectralFD::~ConvFaceTermRHSSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::execute()
{
  CFTRACEBEGIN;
// CF_DEBUG_EXIT;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

  // set the data needed to compute the face terms;
  setFaceTermData();

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

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

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      // if one of the neighbouring cells is parallel updatable or the gradients have to be computed,
      // set face data and reconstruct states
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable() || hasDiffTerm)
      {
        // set the current face and compute the face data in the face term computer
        m_faceTermComputer->setCurrentFace(m_face);
        m_faceTermComputer->computeFaceData();

        // reconstruct the states
        m_faceTermComputer->reconstructFluxPntsStates(m_states);
      }

      // if one of the neighbouring cells is parallel updatable,
      // compute the face term
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // compute the face terms and the wave speed updates
        m_faceTermComputer
            ->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_waveSpeedUpd);

        // update the rhs
        updateRHS();

        // update the wave speeds
        updateWaveSpeed();
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradientFaceTerm();
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
// CF_DEBUG_EXIT;
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::setFaceTermData()
{
  // set the face term data in the face term computer
  m_faceTermComputer->setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::updateRHS()
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
    for (CFuint iRes = 0; iRes < nbrRes; ++iRes, ++resID)
    {
      rhs[resID] += resFactor*m_resUpdates[iSide][iRes];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();
      updateCoeff[solID] += m_waveSpeedUpd[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::computeGradientFaceTerm()
{
  // compute the face term contribution to the gradients
  m_faceTermComputer->computeGradientFaceTerm(m_gradUpdates);

  // add updates to gradients
  addGradFaceTerms();
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::addGradFaceTerms()
{
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_gradUpdates[iSide].size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();

      // update gradients
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradients[solID][iGrad] += m_gradUpdates[iSide][iSol][iGrad];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::setup()
{
  CFAUTOTRACE;

  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_faceTermComputer = getMethodData().getFaceTermComputer();

  // dimensionality and number of equations
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();
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
  m_states.resize(2);
  m_waveSpeedUpd.resize(2);
  m_cells.resize(2);

  m_gradUpdates.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_gradUpdates[iSide].resize(maxNbrSolPnts);
    for (CFuint iSol = 0; iSol < maxNbrSolPnts; ++iSol)
    {
      m_gradUpdates[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_gradUpdates[iSide][iSol][iVar].resize(dim);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFD::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvFaceTermRHSSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
