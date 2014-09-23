#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFV/DiffFaceTermRHSSpectralFV.hh"
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

MethodCommandProvider<DiffFaceTermRHSSpectralFV, SpectralFVMethodData, SpectralFVModule> DiffFaceTermRHSSpectralFVProvider("DiffFaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

DiffFaceTermRHSSpectralFV::DiffFaceTermRHSSpectralFV(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  m_faceTermComputer(CFNULL),
  m_faceBuilder(CFNULL),
  m_face(),
  m_cells(),
  m_cvCVConnSVFace(),
  m_orient(),
  m_states(),
  m_grads(),
  m_resUpdates(),
  m_updateCoefContr(),
  m_nbrEqs()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

DiffFaceTermRHSSpectralFV::~DiffFaceTermRHSSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::configure ( Config::ConfigArgs& args )
{
  SpectralFVMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::execute()
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

      // if one of the neighbouring cells is parallel updatable,
      // compute the face terms
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // set the gradients
        setGradients();

        // set the current face and compute the face data in the face term computer
        m_faceTermComputer->setCurrentFace(m_face);
        m_faceTermComputer->computeFaceData();
        m_faceTermComputer->computeNeighbourCellData();

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
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
// CF_DEBUG_EXIT;
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::setFaceTermData()
{
  // set the face term data in the face term computer
  m_faceTermComputer->setFaceTermData();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get CV-CV connectivity through SV face
  m_cvCVConnSVFace = svLocalData[0]->getExtSVFaceCVConnPerOrient();
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::setGradients()
{
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrStates = m_states[iSide]->size();
    m_grads[iSide].resize(nbrStates);
    for (CFuint iState = 0; iState < nbrStates; ++iState)
    {
      const CFuint stateID = (*m_states[iSide])[iState]->getLocalID();
      m_grads[iSide][iState] = &gradients[stateID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::updateRHS()
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
      rhs[resIDL] += update;
      rhs[resIDR] -= update;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::addUpdateCoeffContributions()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrCVs = m_states[iSide]->size();
    for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
    {
      const CFuint cvID = (*m_states[iSide])[iCV]->getLocalID();
      updateCoeff[cvID] += m_updateCoefContr[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::setup()
{
  CFAUTOTRACE;

  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_faceTermComputer = getMethodData().getFaceTermComputer();

  // number of equations
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
  m_grads.resize(2);
  m_updateCoefContr.resize(2);
  m_cells.resize(2);
}

//////////////////////////////////////////////////////////////////////////////

void DiffFaceTermRHSSpectralFV::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
DiffFaceTermRHSSpectralFV::needsSockets()
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
