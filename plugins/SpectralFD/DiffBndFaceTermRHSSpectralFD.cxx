#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/DiffBndFaceTermRHSSpectralFD.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndFaceTermRHSSpectralFD, SpectralFDMethodData, SpectralFDModule >
  DiffBndFaceTermRHSSpectralFDProvider("DiffBndFaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

DiffBndFaceTermRHSSpectralFD::DiffBndFaceTermRHSSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  m_faceBuilder(CFNULL),
  m_bndFaceTermComputer(CFNULL),
  m_bcStateComputer(CFNULL),
  m_face(),
  m_intCell(),
  m_orient(),
  m_cellStates(),
  m_cellGrads(),
  m_resUpdates(),
  m_updateCoeffContr(),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndFaceTermRHSSpectralFD::~DiffBndFaceTermRHSSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
DiffBndFaceTermRHSSpectralFD::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::executeOnTrs()
{
  CFAUTOTRACE;

  // set BCStateComputer in the boundary face term computer
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // set the data needed to compute the face terms;
  setFaceTermData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCSpectralFD TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get bndFacesStartIdxs from SpectralFDMethodData
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

        // if cell is parallel updatable, compute the boundary face term
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // set the gradients
          setGradients();

          // set the current face and compute the face data in the boundary face term computer
          m_bndFaceTermComputer->setCurrentFace(m_face);
          m_bndFaceTermComputer->computeFaceData();
          m_bndFaceTermComputer->computeNeighbourCellData();

          // reconstruct the states
          m_bndFaceTermComputer->reconstructFluxPntsStates(*m_cellStates);

          // reconstruct the gradients
          m_bndFaceTermComputer->reconstructFluxPntsGradients(m_cellGrads);

          // compute the face terms and the update coefficient contributions
          m_bndFaceTermComputer->computeDiffFaceTermAndUpdateCoefContributions(m_resUpdates,m_updateCoeffContr);

          // update the rhs
          updateRHS();

          // update the contribution to the update coefficient
          addUpdateCoeffContribution();
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
// CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::setFaceTermData()
{
  // set the face term data in the boundary face term computer
  m_bndFaceTermComputer->setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::setGradients()
{
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  const CFuint nbrStates = m_cellStates->size();
  m_cellGrads.resize(nbrStates);
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    m_cellGrads[iState] = &gradients[stateID];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // add updates to residuals
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrRes = m_resUpdates.size();
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes, ++resID)
  {
    rhs[resID] += resFactor*m_resUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::addUpdateCoeffContribution()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
    updateCoeff[solID] += m_updateCoeffContr;
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::setup()
{
  CFAUTOTRACE;

  SpectralFDMethodCom::setup();

  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_bndFaceTermComputer = getMethodData().getBndFaceTermComputer();

  // dimensionality and number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // maximum number of solution points in a cell
  const CFuint maxNbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  m_resUpdates.resize(maxNbrSolPnts*m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndFaceTermRHSSpectralFD::unsetup()
{
  CFAUTOTRACE;

  SpectralFDMethodCom::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD
