#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFV/ConvFaceTermRHSSpectralFV.hh"
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

MethodCommandProvider<ConvFaceTermRHSSpectralFV, SpectralFVMethodData, SpectralFVModule> ConvFaceTermRHSSpectralFVProvider("ConvFaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermRHSSpectralFV::ConvFaceTermRHSSpectralFV(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  m_faceBuilder(CFNULL),
  m_faceTermComputer(CFNULL),
  m_face(),
  m_cells(),
  m_cvCVConnSVFace(),
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

ConvFaceTermRHSSpectralFV::~ConvFaceTermRHSSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::configure ( Config::ConfigArgs& args )
{
  SpectralFVMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::execute()
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
        m_faceTermComputer->reconstructFluxPntsStates(*m_states[LEFT],*m_states[RIGHT]);
        m_faceTermComputer->reconstructFaceAvgState  (*m_states[LEFT],*m_states[RIGHT]);
      }

      // if one of the neighbouring cells is parallel updatable,
      // compute the face term
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // compute the face terms and the wave speed updates
        m_faceTermComputer
                      ->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,
                                                              m_waveSpeedUpd[LEFT ],
                                                              m_waveSpeedUpd[RIGHT]);

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

void ConvFaceTermRHSSpectralFV::setFaceTermData()
{
  // set the face term data in the face term computer
  m_faceTermComputer->setFaceTermData();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get CV-CV connectivity through SV face
  m_cvCVConnSVFace = svLocalData[0]->getExtSVFaceCVConnPerOrient();
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // add fluxes to residuals
  CFuint iRes = 0;
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get first IDs of residuals in these CVs
    CFuint resIDL = m_nbrEqs*(*m_states[LEFT ])[(*m_cvCVConnSVFace)[m_orient][iSubFace][LEFT ]]->getLocalID();
    CFuint resIDR = m_nbrEqs*(*m_states[RIGHT])[(*m_cvCVConnSVFace)[m_orient][iSubFace][RIGHT]]->getLocalID();

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++iRes, ++resIDL, ++resIDR)
    {
      const CFreal update = resFactor*m_resUpdates[iRes];
      rhs[resIDL] -= update;
      rhs[resIDR] += update;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrCVs = m_states[iSide]->size();
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
    {
      const CFuint cvID = (*m_states[iSide])[iCV]->getLocalID();

      updateCoeff[cvID] += m_waveSpeedUpd[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::computeGradientFaceTerm()
{
  // compute the face term contribution to the gradients
  m_faceTermComputer->computeGradientFaceTerm(m_gradUpdates);

  // add updates to gradients
  addGradFaceTerms();
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::addGradFaceTerms()
{
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // update the gradients
  const CFuint nbrSubFaces = (*m_cvCVConnSVFace)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get internal face neighbour CV IDs in states list
    const CFuint cvIDL = (*m_states[LEFT ])[(*m_cvCVConnSVFace)[m_orient][iSubFace][LEFT ]]->getLocalID();
    const CFuint cvIDR = (*m_states[RIGHT])[(*m_cvCVConnSVFace)[m_orient][iSubFace][RIGHT]]->getLocalID();

    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[cvIDL][iGrad] += m_gradUpdates[iSubFace][iGrad];
      gradients[cvIDR][iGrad] -= m_gradUpdates[iSubFace][iGrad];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::setup()
{
  CFAUTOTRACE;

  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_faceTermComputer = getMethodData().getFaceTermComputer();

  // dimensionality and number of equations
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  cf_assert(svLocalData.size() > 0);

  // get SV face - CV connectivity
  SafePtr< vector< vector< CFuint > > > svFaceCVConn = svLocalData[0]->getExtSVFaceCVConn();

  // maximum number of CVs at SV face
  CFuint maxNbrCVsAtSVFace = 0;
  for (CFuint iFace = 0; iFace < svFaceCVConn->size(); ++iFace)
  {
    const CFuint nbrCVsAtSVFace = (*svFaceCVConn)[iFace].size();
    maxNbrCVsAtSVFace = maxNbrCVsAtSVFace > nbrCVsAtSVFace ? maxNbrCVsAtSVFace : nbrCVsAtSVFace;
  }

  // resize variables
  m_resUpdates.resize(maxNbrCVsAtSVFace*m_nbrEqs);
  m_states.resize(2);
  m_waveSpeedUpd.resize(2);
  m_cells.resize(2);

  m_gradUpdates.resize(maxNbrCVsAtSVFace);
  for (CFuint iCV = 0; iCV < maxNbrCVsAtSVFace; ++iCV)
  {
    m_gradUpdates[iCV].resize(m_nbrEqs);
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_gradUpdates[iCV][iVar].resize(dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermRHSSpectralFV::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvFaceTermRHSSpectralFV::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
