#include "Framework/MethodCommandProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/ConvBndFaceTermRHSSpectralFV.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndFaceTermRHSSpectralFV, SpectralFVMethodData, SpectralFVModule >
  ConvBndFaceTermRHSSpectralFVProvider("ConvBndFaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSSpectralFV::ConvBndFaceTermRHSSpectralFV(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  m_faceBuilder(CFNULL),
  m_bndFaceTermComputer(CFNULL),
  m_bcStateComputer(CFNULL),
  m_face(),
  m_intCell(),
  m_svFaceCVConn(),
  m_orient(),
  m_cellStates(),
  m_resUpdates(),
  m_waveSpeedUpd(),
  m_gradUpdates(),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSSpectralFV::~ConvBndFaceTermRHSSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvBndFaceTermRHSSpectralFV::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::configure ( Config::ConfigArgs& args )
{
  SpectralFVMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::executeOnTrs()
{
  CFAUTOTRACE;

  // set BCStateComputer in the boundary face term computer
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // set the data needed to compute the face terms;
  setFaceTermData();

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCSpectralFV TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get bndFacesStartIdxs from SpectralFVMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // set the orientation of the faces
      m_bndFaceTermComputer->setFaceOrientation(m_orient);

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);

        // get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();

        // if cell is parallel updatable or the gradients have to be computed,
        // set face data and reconstruct states
        if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
        {
          // set the current face and compute the face data in the boundary face term computer
          m_bndFaceTermComputer->setCurrentFace(m_face);
          m_bndFaceTermComputer->computeFaceData();

          // reconstruct the states
          m_bndFaceTermComputer->reconstructFluxPntsStates(*m_cellStates);
          m_bndFaceTermComputer->reconstructFaceAvgState  (*m_cellStates);
        }

        // if cell is parallel updatable, compute the boundary face term
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // compute the face terms and the wave speed updates
          m_bndFaceTermComputer->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_waveSpeedUpd);

          // update the rhs
          updateRHS();

          // update the wave speeds in the neighbouring cell
          updateWaveSpeed();
        }

        // if there is a diffusive term, compute the gradients
        if (hasDiffTerm)
        {
          computeGradientFaceTerm();
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
// CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::setFaceTermData()
{
  // set the face term data in the boundary face term computer
  m_bndFaceTermComputer->setFaceTermData();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // boundary face CV connectivity on SV face
  m_svFaceCVConn = svLocalData[0]->getExtSVFaceCVConn();
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // add fluxes to residuals
  CFuint iRes = 0;
  const CFuint nbrSubFaces = (*m_svFaceCVConn)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get first IDs of residual in this CV
    CFuint cvID  = m_nbrEqs*(*m_cellStates)[(*m_svFaceCVConn)[m_orient][iSubFace]]->getLocalID();

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++iRes, ++cvID)
    {
      rhs[cvID] -= resFactor*m_resUpdates[iRes];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbrCVs = m_cellStates->size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    const CFuint cvID = (*m_cellStates)[iCV]->getLocalID();
    updateCoeff[cvID] += m_waveSpeedUpd;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::computeGradientFaceTerm()
{
  // compute the face term contribution to the gradients
  m_bndFaceTermComputer->computeGradientFaceTerm(m_gradUpdates);

  // add updates to gradients
  addGradBCTerms();
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::addGradBCTerms()
{
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // update the gradients
  const CFuint nbrSubFaces = (*m_svFaceCVConn)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // get internal face neighbour CV ID in states list
    const CFuint cvID  = (*m_cellStates)[(*m_svFaceCVConn)[m_orient][iSubFace]]->getLocalID();

    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[cvID][iGrad] += m_gradUpdates[iSubFace][iGrad];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSSpectralFV::setup()
{
  CFAUTOTRACE;

  SpectralFVMethodCom::setup();

  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_bndFaceTermComputer = getMethodData().getBndFaceTermComputer();

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

void ConvBndFaceTermRHSSpectralFV::unsetup()
{
  CFAUTOTRACE;

  SpectralFVMethodCom::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFV

} // namespace COOLFluiD
