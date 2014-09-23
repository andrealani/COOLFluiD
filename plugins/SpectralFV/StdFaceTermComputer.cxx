#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/StdFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    StdFaceTermComputer,SpectralFVMethodData,BaseFaceTermComputer,SpectralFVModule >
  StdFaceTermComputerProvider("StdFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

StdFaceTermComputer::StdFaceTermComputer(const std::string& name) :
  BaseFaceTermComputer(name),
  m_qWheightsSVFaces(),
  m_avgPolyValsSVFaces()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

StdFaceTermComputer::~StdFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void StdFaceTermComputer::setFaceTermData()
{
  // call parent class function
  BaseFaceTermComputer::setFaceTermData();

  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get quadrature wheights on SV faces
  m_qWheightsSVFaces = svLocalData[0]->getExtFaceQWheightsPerOrient();

  // get quadrature node polynomial values on SV faces
  m_svFaceflxPntsRecCoefs = svLocalData[0]->getExtQPntPolyValsPerOrient();

  // get polynomial values averaged over SV face centers
  m_avgPolyValsSVFaces = svLocalData[0]->getSVFaceAvgPolyValsPerOrient();

  // set number of flux points
  const CFuint nbrOrients = m_svFaceflxPntsRecCoefs->size();
  m_nbrFlxPnts.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_nbrFlxPnts[iOrient] = (*m_svFaceflxPntsRecCoefs)[iOrient][LEFT].size();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdFaceTermComputer::reconstructFaceAvgState(const vector< State* >& cellLStates,
                                                  const vector< State* >& cellRStates)
{
  m_statesReconstr->reconstructState(cellLStates,*m_faceAvgSol[LEFT ],
                                     (*m_avgPolyValsSVFaces)[m_orient][LEFT ],
                                     cellLStates.size());
  m_statesReconstr->reconstructState(cellRStates,*m_faceAvgSol[RIGHT],
                                     (*m_avgPolyValsSVFaces)[m_orient][RIGHT],
                                     cellRStates.size());
}

//////////////////////////////////////////////////////////////////////////////

void StdFaceTermComputer::computeFaceFluxIntegralFromFlxPntFluxes(RealVector& resUpdates)
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
    const CFuint resIDStore = resID;

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

void StdFaceTermComputer::computeGradFaceTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
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

void StdFaceTermComputer::setup()
{
  CFAUTOTRACE;

  BaseFaceTermComputer::setup();

  // dimensionality and number of variables
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  cf_assert(svLocalData.size() > 0);

  // number of external quadrature nodes
  SafePtr< vector< vector< CFreal > > > qWheightsSVFaces = svLocalData[0]->getExtQWheightsPerOrient();
  cf_assert(qWheightsSVFaces->size() > 0);
  const CFuint nbrQNodes = (*qWheightsSVFaces)[0].size();

  // resize variable for unit normals in quadrature nodes
  m_unitNormalFlxPnt.resize(nbrQNodes);
  for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
  {
    m_unitNormalFlxPnt[iQNode].resize(dim);
  }

  // create left and right states and extra variables
  m_flxPntSol      .resize(2);
  m_flxPntExtraVars.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_flxPntSol      [iSide].resize(nbrQNodes);
    m_flxPntExtraVars[iSide].resize(nbrQNodes);
    for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
    {
      m_flxPntSol      [iSide][iQNode] = new State();
      m_flxPntExtraVars[iSide][iQNode] = new RealVector(m_nbrExtraVars);
    }
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(nbrQNodes);

  // set m_allSol and m_allExtraVars pointers to left and right states and extra variables
  m_allSol      .resize(2*(nbrQNodes+1));
  m_allExtraVars.resize(2*(nbrQNodes+1));
  CFuint iState = 0;
  for (CFuint iQNode = 0; iQNode < nbrQNodes; ++iQNode)
  {
    m_allSol      [iState] = m_flxPntSol      [LEFT ][iQNode];
    m_allExtraVars[iState] = m_flxPntExtraVars[LEFT ][iQNode]; ++iState;
    m_allSol      [iState] = m_flxPntSol      [RIGHT][iQNode];
    m_allExtraVars[iState] = m_flxPntExtraVars[RIGHT][iQNode]; ++iState;
  }
  m_allSol      [iState] = m_faceAvgSol      [LEFT ];
  m_allExtraVars[iState] = m_faceAvgExtraVars[LEFT ]; ++iState;
  m_allSol      [iState] = m_faceAvgSol      [RIGHT];
  m_allExtraVars[iState] = m_faceAvgExtraVars[RIGHT];

  // setup variables for gradient computation
  if (getMethodData().hasDiffTerm())
  {
    // create gradient vars for quadrature points
    m_flxPntGradVars.resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_flxPntGradVars[iSide].resize(m_nbrEqs);
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntGradVars[iSide][iVar] = new RealVector(nbrQNodes);
      }
    }

    // set RealVector pointers to states
    m_flxPntRVSol.resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_flxPntRVSol[iSide].resize(nbrQNodes);
      for (CFuint iQPnt = 0; iQPnt < nbrQNodes; ++iQPnt)
      {
        m_flxPntRVSol[iSide][iQPnt] = m_flxPntSol[iSide][iQPnt];
      }
    }

    // create gradients for quadrature points
    m_flxPntGrads.resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_flxPntGrads[iSide].resize(nbrQNodes);
      for (CFuint iQPnt = 0; iQPnt < nbrQNodes; ++iQPnt)
      {
        m_flxPntGrads[iSide][iQPnt].resize(m_nbrEqs);
        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
        {
          m_flxPntGrads[iSide][iQPnt][iGrad] = new RealVector(dim);
        }
      }
    }

    // create pointers to gradients for quadrature points
    m_flxPntGradPtrs.resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_flxPntGradPtrs[iSide].resize(nbrQNodes);
      for (CFuint iQPnt = 0; iQPnt < nbrQNodes; ++iQPnt)
      {
        m_flxPntGradPtrs[iSide][iQPnt] = &m_flxPntGrads[iSide][iQPnt];
      }
    }
  }

  // set maximum number of quadrature points that will be passed at one time to the physical model
  const CFuint maxNbrNodes = 2*(nbrQNodes+1);
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

void StdFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSide = 0; iSide < m_flxPntSol.size(); ++iSide)
  {
    for (CFuint iQNode = 0; iQNode < m_flxPntSol[iSide].size(); ++iQNode)
    {
      deletePtr(m_flxPntSol      [iSide][iQNode]);
      deletePtr(m_flxPntExtraVars[iSide][iQNode]);
    }
    m_flxPntSol      [iSide].resize(0);
    m_flxPntExtraVars[iSide].resize(0);
  }
  m_flxPntSol      .resize(0);
  m_flxPntExtraVars.resize(0);

  for (CFuint iSide = 0; iSide < m_flxPntGradVars.size(); ++iSide)
  {
    for (CFuint iState = 0; iState < m_flxPntGradVars[iSide].size(); ++iState)
    {
      deletePtr(m_flxPntGradVars[iSide][iState]);
    }
    m_flxPntGradVars[iSide].resize(0);
  }
  m_flxPntGradVars.resize(0);

  for (CFuint iSide = 0; iSide < m_flxPntGrads.size(); ++iSide)
  {
    for (CFuint iQPnt = 0; iQPnt < m_flxPntGrads[iSide].size(); ++iQPnt)
    {
      for (CFuint iGrad = 0; iGrad < m_flxPntGrads[iSide][iQPnt].size(); ++iGrad)
      {
        deletePtr(m_flxPntGrads[iSide][iQPnt][iGrad]);
      }
      m_flxPntGrads[iSide][iQPnt].resize(0);
    }
    m_flxPntGrads[iSide].resize(0);
  }
  m_flxPntGrads.resize(0);

  BaseFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
