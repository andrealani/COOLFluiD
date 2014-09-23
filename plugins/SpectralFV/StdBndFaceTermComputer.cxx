#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/StdBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    StdBndFaceTermComputer,SpectralFVMethodData,BaseBndFaceTermComputer,SpectralFVModule >
  StdBndFaceTermComputerProvider("StdBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

StdBndFaceTermComputer::StdBndFaceTermComputer(const std::string& name) :
  BaseBndFaceTermComputer(name),
  m_qWheightsSVFaces(),
  m_avgPolyValsSVFaces()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

StdBndFaceTermComputer::~StdBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void StdBndFaceTermComputer::setFaceTermData()
{
  // call parent class function
  BaseBndFaceTermComputer::setFaceTermData();

  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get quadrature wheights on SV faces
  m_qWheightsSVFaces = svLocalData[0]->getExtFaceQuadWheights();

  // get quadrature node wheight coordinates on SV faces
  m_flxPntWheightCoordsSVFaces = svLocalData[0]->getExtQPntWheightCoords();

  // get quadrature node polynomial values on SV faces
  m_svFaceflxPntsRecCoefs = svLocalData[0]->getExtQPntPolyVals();

  // get polynomial values averaged over SV face centers
  m_avgPolyValsSVFaces = svLocalData[0]->getSVFaceAvgPolyVals();

  // set number of flux points
  const CFuint nbrOrients = m_svFaceflxPntsRecCoefs->size();
  m_nbrFlxPnts.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_nbrFlxPnts[iOrient] = (*m_svFaceflxPntsRecCoefs)[iOrient].size();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdBndFaceTermComputer::reconstructFaceAvgState(const vector< State* >& cellIntStates)
{
  m_statesReconstr->reconstructState(cellIntStates,*m_faceAvgSolInt,
                                     (*m_avgPolyValsSVFaces)[m_orient],
                                     cellIntStates.size());
}

//////////////////////////////////////////////////////////////////////////////

void StdBndFaceTermComputer::computeFaceFluxIntegralFromFlxPntFluxes(RealVector& resUpdates)
{
  // set residual updates to zero
  resUpdates = 0.0;

  // compute face term
  CFuint stateID = 0;
  CFuint resID = 0;
  const CFuint nbrSubFaces = (*m_qWheightsSVFaces)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // store resID
    CFuint resIDStore = resID;

    // Calculate the Riemann flux integrated over the face
    const CFuint nbrExtFaceQPnts = (*m_qWheightsSVFaces)[m_orient][iSubFace].size();
    for (CFuint iQPnt = 0; iQPnt < nbrExtFaceQPnts; ++iQPnt, ++stateID)
    {
      resID = resIDStore;
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resID)
      {
        // add local Riemann flux to total flux through face
        resUpdates[resID] += (*m_qWheightsSVFaces)[m_orient][iSubFace][iQPnt]*m_flxPntRiemannFlx[stateID][iVar];
      }
    }
  }
  resUpdates *= m_surf;
}

//////////////////////////////////////////////////////////////////////////////

void StdBndFaceTermComputer::computeGradFaceTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
{
  CFuint stateID = 0;
  const CFuint nbrSubFaces = (*m_qWheightsSVFaces)[m_orient].size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    // set m_gradUpdates[iSubFace] to zero
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradUpdates[iSubFace][iGrad] = 0.0;
    }

    // Calculate the Riemann flux integrated over the face
    const CFuint nbrExtFaceQPnts = (*m_qWheightsSVFaces)[m_orient][iSubFace].size();
    for (CFuint iQPnt = 0; iQPnt < nbrExtFaceQPnts; ++iQPnt, ++stateID)
    {
      const CFreal qWheight = (*m_qWheightsSVFaces)[m_orient][iSubFace][iQPnt];
      RealVector& normal = m_unitNormalFlxPnt[stateID];
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        // add contribution to gradient term
        gradUpdates[iSubFace][iGrad] += qWheight*m_flxPntRiemannFlx[stateID][iGrad]*normal;
      }
    }

    // multiply updates with face surface
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradUpdates[iSubFace][iGrad] *= m_surf;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  BaseBndFaceTermComputer::setup();

  // dimensionality and number of variables
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  cf_assert(svLocalData.size() > 0);

  // number of external quadrature nodes
  SafePtr< vector< vector< CFreal > > > qWheightsSVFaces = svLocalData[0]->getExtQWheightsPerOrient();
  cf_assert(qWheightsSVFaces->size() > 0);
  const CFuint nbrQNodes = (*qWheightsSVFaces)[0].size();

  // resize variables
  m_unitNormalFlxPnt.resize(nbrQNodes);
  for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
  {
    m_unitNormalFlxPnt[iQNode].resize(dim);
  }

  m_flxPntCoords.resize(nbrQNodes);
  for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
  {
    m_flxPntCoords[iQNode].resize(dim);
  }

  // create internal and ghost states and extra variables
  m_flxPntIntSol        .resize(nbrQNodes);
  m_flxPntGhostSol      .resize(nbrQNodes);
  m_flxPntIntExtraVars  .resize(nbrQNodes);
  m_flxPntGhostExtraVars.resize(nbrQNodes);
  for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
  {
    m_flxPntIntSol        [iQNode] = new State();
    m_flxPntGhostSol      [iQNode] = new State();
    m_flxPntIntExtraVars  [iQNode] = new RealVector(m_nbrExtraVars);
    m_flxPntGhostExtraVars[iQNode] = new RealVector(m_nbrExtraVars);
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(nbrQNodes);

  // resize m_backupGhostStates
  m_backupGhostStates.resize(nbrQNodes);
  for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
  {
    m_backupGhostStates[iQNode].resize(m_nbrEqs);
  }

  // set m_allSol and m_allExtraVars pointers
  m_allSol      .resize(2*nbrQNodes+1);
  m_allExtraVars.resize(2*nbrQNodes+1);
  CFuint iState = 0;
  for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
  {
    m_allSol      [iState] = m_flxPntIntSol        [iQNode];
    m_allExtraVars[iState] = m_flxPntIntExtraVars  [iQNode]; ++iState;
    m_allSol      [iState] = m_flxPntGhostSol      [iQNode];
    m_allExtraVars[iState] = m_flxPntGhostExtraVars[iQNode]; ++iState;
  }
  m_allSol      [iState] = m_faceAvgSolInt     ;
  m_allExtraVars[iState] = m_faceAvgExtraVarInt;

  // setup variables for gradient computation
  if (getMethodData().hasDiffTerm())
  {
    // create gradient vars for quadrature points
    m_flxPntIntGradVars  .resize(m_nbrEqs);
    m_flxPntGhostGradVars.resize(m_nbrEqs);
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_flxPntIntGradVars  [iVar] = new RealVector(nbrQNodes);
      m_flxPntGhostGradVars[iVar] = new RealVector(nbrQNodes);
    }

    // set RealVector pointers to states
    m_flxPntIntRVSol  .resize(nbrQNodes);
    m_flxPntGhostRVSol.resize(nbrQNodes);
    for (CFuint iQPnt = 0; iQPnt < nbrQNodes; ++iQPnt)
    {
      m_flxPntIntRVSol  [iQPnt] = m_flxPntIntSol  [iQPnt];
      m_flxPntGhostRVSol[iQPnt] = m_flxPntGhostSol[iQPnt];
    }

    // create gradients for quadrature points
    m_flxPntIntGrads  .resize(nbrQNodes);
    m_flxPntGhostGrads.resize(nbrQNodes);
    for (CFuint iQPnt = 0; iQPnt < nbrQNodes; ++iQPnt)
    {
      m_flxPntIntGrads  [iQPnt].resize(m_nbrEqs);
      m_flxPntGhostGrads[iQPnt].resize(m_nbrEqs);
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        m_flxPntIntGrads  [iQPnt][iGrad] = new RealVector(dim);
        m_flxPntGhostGrads[iQPnt][iGrad] = new RealVector(dim);
      }
    }

    // create pointers to gradients for quadrature points
    m_flxPntIntGradPtrs  .resize(nbrQNodes);
    m_flxPntGhostGradPtrs.resize(nbrQNodes);
    for (CFuint iQPnt = 0; iQPnt < nbrQNodes; ++iQPnt)
    {
      m_flxPntIntGradPtrs  [iQPnt] = &m_flxPntIntGrads  [iQPnt];
      m_flxPntGhostGradPtrs[iQPnt] = &m_flxPntGhostGrads[iQPnt];
    }
  }

  // set maximum number of quadrature points that will be passed at one time to the physical model
  const CFuint maxNbrNodes = 2*nbrQNodes+1;
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
    setMaxNbrStatesData(prevMaxNbrStatesData > maxNbrNodes ? prevMaxNbrStatesData : maxNbrNodes);

  // set maximum number of points in which the Riemann flux has to be evaluated at the same time
  const CFuint prevMaxNbrRFluxPnts = getMethodData().getMaxNbrRFluxPnts();
  getMethodData().
    setMaxNbrRFluxPnts
    (
      prevMaxNbrRFluxPnts > nbrQNodes ? prevMaxNbrRFluxPnts : nbrQNodes
    );
}

//////////////////////////////////////////////////////////////////////////////

void StdBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iQNode = 0; iQNode < m_flxPntIntSol.size(); ++iQNode)
  {
    deletePtr(m_flxPntIntSol        [iQNode]);
    deletePtr(m_flxPntGhostSol      [iQNode]);
    deletePtr(m_flxPntIntExtraVars  [iQNode]);
    deletePtr(m_flxPntGhostExtraVars[iQNode]);
  }
  m_flxPntIntSol        .resize(0);
  m_flxPntGhostSol      .resize(0);
  m_flxPntIntExtraVars  .resize(0);
  m_flxPntGhostExtraVars.resize(0);

  for (CFuint iState = 0; iState < m_flxPntIntGradVars.size(); ++iState)
  {
    deletePtr(m_flxPntIntGradVars  [iState]);
    deletePtr(m_flxPntGhostGradVars[iState]);
  }
  m_flxPntIntGradVars  .resize(0);
  m_flxPntGhostGradVars.resize(0);

  for (CFuint iQPnt = 0; iQPnt < m_flxPntIntGrads.size(); ++iQPnt)
  {
    for (CFuint iGrad = 0; iGrad < m_flxPntIntGrads[iQPnt].size(); ++iGrad)
    {
      deletePtr(m_flxPntIntGrads  [iQPnt][iGrad]);
      deletePtr(m_flxPntGhostGrads[iQPnt][iGrad]);
    }
    m_flxPntIntGrads  [iQPnt].resize(0);
    m_flxPntGhostGrads[iQPnt].resize(0);
  }
  m_flxPntIntGrads  .resize(0);
  m_flxPntGhostGrads.resize(0);

  BaseBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
