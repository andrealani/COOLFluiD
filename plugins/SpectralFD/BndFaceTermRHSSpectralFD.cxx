#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/BndFaceTermRHSSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BndFaceTermRHSSpectralFD, SpectralFDMethodData, SpectralFDModule >
  BndFaceTermRHSSpectralFDProvider("BndFaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

BndFaceTermRHSSpectralFD::BndFaceTermRHSSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  m_faceBuilder(CFNULL),
  m_bndFaceTermComputer(CFNULL),
  m_bcStateComputer(CFNULL),
  m_face(),
  m_intCell(),
  m_orient(),
  m_cellStates(),
  m_resUpdates(),
  m_diffResUpdates(),
  m_updateCoeffContr(),
  m_diffUpdateCoeffContr(),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

BndFaceTermRHSSpectralFD::~BndFaceTermRHSSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
BndFaceTermRHSSpectralFD::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSSpectralFD::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSSpectralFD::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current TRS
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

  // SET THE DATA NEEDED TO COMPUTE THE FACE TERMS
  setFaceTermComputerData();

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // SET THE ORIENTATION OF THE FACES
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
          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceTermData();

          // COMPUTE CONVECTIVE BOUNDARY FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
          m_bndFaceTermComputer->
              computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_updateCoeffContr);

          // COMPUTE DIFFUSIVE BOUNDARY FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
          m_bndFaceTermComputer->
              computeDiffFaceTermAndUpdateCoefContributions(m_diffResUpdates,m_diffUpdateCoeffContr);

          // ADD TOTAL UPDATE TO RESIDUAL AND UPDATE COEFFICIENTS
          m_resUpdates += m_diffResUpdates;
          addUpdatesToResidual();
          m_updateCoeffContr += m_diffUpdateCoeffContr;
          addUpdateCoeffContribution();
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSSpectralFD::setFaceTermComputerData()
{
  // set BCStateComputer in the boundary face term computer
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // set the face term data in the boundary face term computer
  m_bndFaceTermComputer->setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSSpectralFD::computeFaceTermData()
{
  // set the face in the boundary face term computer
  m_bndFaceTermComputer->setCurrentFace(m_face);

  // compute the face data in the boundary face term computer
  m_bndFaceTermComputer->computeFaceData();

  // compute the neighbouring cell data in the boundary face term computer
  m_bndFaceTermComputer->computeNeighbourCellData();

  // reconstruct the states in the face flux points
  m_bndFaceTermComputer->reconstructFluxPntsStates(*m_cellStates);

  // compute the solution polynomial gradients in the face flux points
  m_bndFaceTermComputer->reconstructFluxPntsSolPolyGrads(*m_cellStates);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSSpectralFD::addUpdatesToResidual()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // add updates to residuals
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrRes = m_resUpdates.size();
  /// @warning this size is for the highest polynomial degree, this should probably be changed for p-multigrid
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes, ++resID)
  {
    rhs[resID] += resFactor*m_resUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSSpectralFD::addUpdateCoeffContribution()
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

void BndFaceTermRHSSpectralFD::setup()
{
  CFAUTOTRACE;

  SpectralFDMethodCom::setup();

  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_bndFaceTermComputer = getMethodData().getBndFaceTermComputer().d_castTo<CompactBndFaceTermComputer>();

  // dimensionality and number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // maximum number of solution points in a cell
  const CFuint maxNbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  m_resUpdates    .resize(maxNbrSolPnts*m_nbrEqs);
  m_diffResUpdates.resize(maxNbrSolPnts*m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void BndFaceTermRHSSpectralFD::unsetup()
{
  CFAUTOTRACE;

  SpectralFDMethodCom::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD
