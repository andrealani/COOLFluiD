#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMethod/ConvBndFaceTermRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndFaceTermRHSFluxReconstruction, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndFaceTermRHSFluxReconstructionProvider("ConvBndFaceTermRHS");

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSFluxReconstruction::ConvBndFaceTermRHSFluxReconstruction(const std::string& name) :
  FluxReconstructionSolverCom(name),
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
  m_resUpdates(),
  m_waveSpeedUpd(),
  m_gradUpdates(),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndFaceTermRHSFluxReconstruction::~ConvBndFaceTermRHSFluxReconstruction()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvBndFaceTermRHSFluxReconstruction::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSFluxReconstruction::executeOnTrs()
{
  CFAUTOTRACE;

  // set BCStateComputer in the boundary face term computer
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // set the data needed to compute the face terms;
  setFaceTermData();

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = false;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCFluxReconstruction TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get bndFacesStartIdxs from FluxReconstructionMethodData
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

void ConvBndFaceTermRHSFluxReconstruction::setFaceTermData()
{
  // set the face term data in the boundary face term computer
  m_bndFaceTermComputer->setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSFluxReconstruction::updateRHS()
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

void ConvBndFaceTermRHSFluxReconstruction::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
    updateCoeff[solID] += m_waveSpeedUpd;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSFluxReconstruction::computeGradientFaceTerm()
{
  // compute the face term contribution to the gradients
  m_bndFaceTermComputer->computeGradientFaceTerm(m_gradUpdates);

  // add updates to gradients
  addGradBCTerms();
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSFluxReconstruction::addGradBCTerms()
{
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // update the gradients
  const CFuint nbrSolPnts = m_gradUpdates.size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
      // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

      // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[iSol][iGrad];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSFluxReconstruction::setup()
{
  CFAUTOTRACE;

  FluxReconstructionSolverCom::setup();

  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get and setup the face term computer
  m_bndFaceTermComputer = getMethodData().getBndFaceTermComputer();

  // dimensionality and number of equations
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& sdLocalData = getMethodData().getFRLocalData();
  cf_assert(sdLocalData.size() > 0);

  // maximum number of solution points in a cell
  const CFuint maxNbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  m_resUpdates.resize(maxNbrSolPnts*m_nbrEqs);

  m_gradUpdates.resize(maxNbrSolPnts);
  for (CFuint iSol = 0; iSol < maxNbrSolPnts; ++iSol)
  {
    m_gradUpdates[iSol].resize(m_nbrEqs);
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_gradUpdates[iSol][iVar].resize(dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndFaceTermRHSFluxReconstruction::unsetup()
{
  CFAUTOTRACE;

  FluxReconstructionSolverCom::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
