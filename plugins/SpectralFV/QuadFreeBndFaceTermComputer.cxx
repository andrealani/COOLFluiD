#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/QuadFreeBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    QuadFreeBndFaceTermComputer,SpectralFVMethodData,BaseBndFaceTermComputer,SpectralFVModule >
  QuadFreeBndFaceTermComputerProvider("QuadFreeBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

QuadFreeBndFaceTermComputer::QuadFreeBndFaceTermComputer(const std::string& name) :
  BaseBndFaceTermComputer(name),
  m_cvExtFaceFluxCoef(),
  m_avgSolInSVFaceCoef()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

QuadFreeBndFaceTermComputer::~QuadFreeBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeBndFaceTermComputer::setFaceTermData()
{
  // call parent class function
  BaseBndFaceTermComputer::setFaceTermData();

  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get coefficients for the reconstruction of the solution in the face flux points
  m_svFaceflxPntsRecCoefs = svLocalData[0]->getSolInFaceFluxPntCoef();

  // get quadrature node coordinates on SV faces
  m_flxPntWheightCoordsSVFaces = svLocalData[0]->getFaceFluxPolyNodeWheightCoord();

  // get coefficients for the computation of the flux through an external CV face
  m_cvExtFaceFluxCoef = svLocalData[0]->getCVExtFaceFluxCoef();

  // get coefficients for the computation of the average solution over an external CV face
  m_avgSolInSVFaceCoef = svLocalData[0]->getAvgSolInSVFaceCoef();

  // set number of flux points
  const CFuint nbrOrients = m_svFaceflxPntsRecCoefs->size();
  m_nbrFlxPnts.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_nbrFlxPnts[iOrient] = (*m_svFaceflxPntsRecCoefs)[iOrient].size();
  }
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeBndFaceTermComputer::reconstructFaceAvgState(const vector< State* >& cellIntStates)
{
  m_statesReconstr->reconstructState(m_flxPntIntSol,*m_faceAvgSolInt,
                                     *m_avgSolInSVFaceCoef,
                                     m_flxPntIntSol.size());
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeBndFaceTermComputer::computeFaceFluxIntegralFromFlxPntFluxes(RealVector& resUpdates)
{
  // set residual updates to zero
  resUpdates = 0.0;

  // compute face term
  CFuint resID = 0;
  const CFuint nbrSubFaces = m_cvExtFaceFluxCoef->size();
  const CFuint nbrFaceFlxPnts = m_flxPntRiemannFlx.size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubFaces; ++iSubFace)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resID)
    {
      for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
      {
        // constribution to flux through cv faces
        resUpdates[resID] += (*m_cvExtFaceFluxCoef)[iSubFace][iFlx]*m_flxPntRiemannFlx[iFlx][iVar];
      }
    }
  }
  // multiply with face surface
  resUpdates *= m_surf;
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeBndFaceTermComputer::computeGradFaceTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
{
  const CFuint nbrSubfaces = m_cvExtFaceFluxCoef->size();
  const CFuint nbrFaceFlxPnts = m_flxPntIntRVSol.size();
  for (CFuint iSubFace = 0; iSubFace < nbrSubfaces; ++iSubFace)
  {
    // set m_gradUpdates[iSubFace] to zero
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradUpdates[iSubFace][iGrad] = 0.0;
    }

    // Calculate the Riemann flux integrated over the face
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal flxPntWheight = (*m_cvExtFaceFluxCoef)[iSubFace][iFlx];
      RealVector& normal = m_unitNormalFlxPnt[iFlx];
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        // add contribution to gradient term
        gradUpdates[iSubFace][iGrad] += flxPntWheight*m_flxPntRiemannFlx[iFlx][iGrad]*normal;
      }
    }

    // multiply with face surface
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradUpdates[iSubFace][iGrad] *= m_surf;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  BaseBndFaceTermComputer::setup();

  // dimensionality and number of variables
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();

  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get maximum number of flux points at a SV face
  const CFuint maxNbrFlxPntsFaces = svLocalData[0]->getNbrOfFaceFlxPnts();

  // resize variable for unit normals in flux points
  m_unitNormalFlxPnt.resize(maxNbrFlxPntsFaces);
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
  {
    m_unitNormalFlxPnt[iFlx].resize(dim);
  }

  m_flxPntCoords.resize(maxNbrFlxPntsFaces);
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(dim);
  }

  // create flux point states and extra variables
  m_flxPntIntSol        .resize(maxNbrFlxPntsFaces);
  m_flxPntGhostSol      .resize(maxNbrFlxPntsFaces);
  m_flxPntIntExtraVars  .resize(maxNbrFlxPntsFaces);
  m_flxPntGhostExtraVars.resize(maxNbrFlxPntsFaces);
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
  {
    m_flxPntIntSol        [iFlx] = new State();
    m_flxPntGhostSol      [iFlx] = new State();
    m_flxPntIntExtraVars  [iFlx] = new RealVector(m_nbrExtraVars);
    m_flxPntGhostExtraVars[iFlx] = new RealVector(m_nbrExtraVars);
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(maxNbrFlxPntsFaces);

    // resize m_backupGhostStates
  m_backupGhostStates.resize(maxNbrFlxPntsFaces);
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
  {
    m_backupGhostStates[iFlx].resize(m_nbrEqs);
  }

  // set m_allSol pointers to left and right states
  const CFuint maxNbrFlxPntsFacesX2P1 = 2*maxNbrFlxPntsFaces+1;
  m_allSol      .resize(maxNbrFlxPntsFacesX2P1);
  m_allExtraVars.resize(maxNbrFlxPntsFacesX2P1);
  CFuint iState = 0;
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
  {
    m_allSol      [iState] = m_flxPntIntSol        [iFlx];
    m_allExtraVars[iState] = m_flxPntIntExtraVars  [iFlx]; ++iState;
    m_allSol      [iState] = m_flxPntGhostSol      [iFlx];
    m_allExtraVars[iState] = m_flxPntGhostExtraVars[iFlx]; ++iState;
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
      m_flxPntIntGradVars  [iVar] = new RealVector(maxNbrFlxPntsFaces);
      m_flxPntGhostGradVars[iVar] = new RealVector(maxNbrFlxPntsFaces);
    }

    // set RealVector pointers to states
    m_flxPntIntRVSol  .resize(maxNbrFlxPntsFaces);
    m_flxPntGhostRVSol.resize(maxNbrFlxPntsFaces);
    for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
    {
      m_flxPntIntRVSol  [iFlx] = m_flxPntIntSol  [iFlx];
      m_flxPntGhostRVSol[iFlx] = m_flxPntGhostSol[iFlx];
    }

    // create gradients for flux points
    m_flxPntIntGrads  .resize(maxNbrFlxPntsFaces);
    m_flxPntGhostGrads.resize(maxNbrFlxPntsFaces);
    for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
    {
      m_flxPntIntGrads  [iFlx].resize(m_nbrEqs);
      m_flxPntGhostGrads[iFlx].resize(m_nbrEqs);
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        m_flxPntIntGrads  [iFlx][iGrad] = new RealVector(dim);
        m_flxPntGhostGrads[iFlx][iGrad] = new RealVector(dim);
      }
    }

    // create pointers to gradients for flux points
    m_flxPntIntGradPtrs  .resize(maxNbrFlxPntsFaces);
    m_flxPntGhostGradPtrs.resize(maxNbrFlxPntsFaces);
    for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
    {
      m_flxPntIntGradPtrs  [iFlx] = &m_flxPntIntGrads  [iFlx];
      m_flxPntGhostGradPtrs[iFlx] = &m_flxPntGhostGrads[iFlx];
    }
  }

  // set maximum number of flux points that will be passed at one time to the physical model
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
    setMaxNbrStatesData
    (
      prevMaxNbrStatesData > maxNbrFlxPntsFacesX2P1 ? prevMaxNbrStatesData : maxNbrFlxPntsFacesX2P1
    );

  // set maximum number of points in which the Riemann flux has to be evaluated at the same time
  const CFuint prevMaxNbrRFluxPnts = getMethodData().getMaxNbrRFluxPnts();
  getMethodData().
    setMaxNbrRFluxPnts
    (
      prevMaxNbrRFluxPnts > maxNbrFlxPntsFaces ? prevMaxNbrRFluxPnts : maxNbrFlxPntsFaces
    );
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iFlx = 0; iFlx < m_flxPntIntSol.size(); ++iFlx)
  {
    deletePtr(m_flxPntIntSol        [iFlx]);
    deletePtr(m_flxPntGhostSol      [iFlx]);
    deletePtr(m_flxPntIntExtraVars  [iFlx]);
    deletePtr(m_flxPntGhostExtraVars[iFlx]);
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

  for (CFuint iState = 0; iState < m_flxPntIntGrads.size(); ++iState)
  {
    for (CFuint iGrad = 0; iGrad < m_flxPntIntGrads[iState].size(); ++iGrad)
    {
      deletePtr(m_flxPntIntGrads  [iState][iGrad]);
      deletePtr(m_flxPntGhostGrads[iState][iGrad]);
    }
    m_flxPntIntGrads  [iState].resize(0);
    m_flxPntGhostGrads[iState].resize(0);
  }
  m_flxPntIntGrads  .resize(0);
  m_flxPntGhostGrads.resize(0);

  BaseBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
