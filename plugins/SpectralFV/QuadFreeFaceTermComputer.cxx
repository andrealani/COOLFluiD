#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/QuadFreeFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    QuadFreeFaceTermComputer,SpectralFVMethodData,BaseFaceTermComputer,SpectralFVModule >
  QuadFreeFaceTermComputerProvider("QuadFreeFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

QuadFreeFaceTermComputer::QuadFreeFaceTermComputer(const std::string& name) :
  BaseFaceTermComputer(name),
  m_cvExtFaceFluxCoef(),
  m_avgSolInSVFaceCoef()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

QuadFreeFaceTermComputer::~QuadFreeFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeFaceTermComputer::setFaceTermData()
{
  // call parent class function
  BaseFaceTermComputer::setFaceTermData();

  // get the SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get coefficients for the reconstruction of the solution in the face flux points
  m_svFaceflxPntsRecCoefs = svLocalData[0]->getSolInFaceFluxPntCoefPerOrient();

  // get coefficients for the computation of the flux through an external CV face
  m_cvExtFaceFluxCoef = svLocalData[0]->getCVExtFaceFluxCoef();

  // get coefficients for the computation of the average solution over an external CV face
  m_avgSolInSVFaceCoef = svLocalData[0]->getAvgSolInSVFaceCoef();

  // set number of flux points
  const CFuint nbrOrients = m_svFaceflxPntsRecCoefs->size();
  m_nbrFlxPnts.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_nbrFlxPnts[iOrient] = (*m_svFaceflxPntsRecCoefs)[iOrient][LEFT].size();
  }
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeFaceTermComputer::reconstructFaceAvgState(const vector< State* >& cellLStates,
                                                       const vector< State* >& cellRStates)
{
  m_statesReconstr->reconstructState(m_flxPntSol[LEFT ],*m_faceAvgSol[LEFT ],
                                     *m_avgSolInSVFaceCoef,
                                     m_flxPntSol[LEFT ].size());
  m_statesReconstr->reconstructState(m_flxPntSol[RIGHT],*m_faceAvgSol[RIGHT],
                                     *m_avgSolInSVFaceCoef,
                                     m_flxPntSol[RIGHT].size());
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeFaceTermComputer::computeFaceFluxIntegralFromFlxPntFluxes(RealVector& resUpdates)
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

void QuadFreeFaceTermComputer::computeGradFaceTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
{
  const CFuint nbrSubfaces = m_cvExtFaceFluxCoef->size();
  const CFuint nbrFaceFlxPnts = m_flxPntRVSol[LEFT].size();
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

void QuadFreeFaceTermComputer::setup()
{
  CFAUTOTRACE;

  BaseFaceTermComputer::setup();

  // dimensionality and number of variables
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

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

  // create flux point states
  m_flxPntSol      .resize(2);
  m_flxPntExtraVars.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_flxPntSol      [iSide].resize(maxNbrFlxPntsFaces);
    m_flxPntExtraVars[iSide].resize(maxNbrFlxPntsFaces);
    for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
    {
      m_flxPntSol      [iSide][iFlx] = new State();
      m_flxPntExtraVars[iSide][iFlx] = new RealVector(m_nbrExtraVars);
    }
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(maxNbrFlxPntsFaces);

  // set m_allSol and m_allExtraVars pointers to left and right states and extra variables
  const CFuint maxNbrFlxPntsFacesX2P2 = 2*maxNbrFlxPntsFaces+2;
  m_allSol      .resize(maxNbrFlxPntsFacesX2P2);
  m_allExtraVars.resize(maxNbrFlxPntsFacesX2P2);
  CFuint iState = 0;
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
  {
    m_allSol      [iState] = m_flxPntSol      [LEFT ][iFlx];
    m_allExtraVars[iState] = m_flxPntExtraVars[LEFT ][iFlx]; ++iState;
    m_allSol      [iState] = m_flxPntSol      [RIGHT][iFlx];
    m_allExtraVars[iState] = m_flxPntExtraVars[RIGHT][iFlx]; ++iState;
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
        m_flxPntGradVars[iSide][iVar] = new RealVector(maxNbrFlxPntsFaces);
      }
    }

    // set RealVector pointers to states
    m_flxPntRVSol.resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_flxPntRVSol[iSide].resize(maxNbrFlxPntsFaces);
      for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
      {
        m_flxPntRVSol[iSide][iFlx] = m_flxPntSol[iSide][iFlx];
      }
    }

    // create gradients for flux points
    m_flxPntGrads.resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_flxPntGrads[iSide].resize(maxNbrFlxPntsFaces);
      for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
      {
        m_flxPntGrads[iSide][iFlx].resize(m_nbrEqs);
        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
        {
          m_flxPntGrads[iSide][iFlx][iGrad] = new RealVector(dim);
        }
      }
    }

    // create pointers to gradients for flux points
    m_flxPntGradPtrs.resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_flxPntGradPtrs[iSide].resize(maxNbrFlxPntsFaces);
      for (CFuint iFlx = 0; iFlx < maxNbrFlxPntsFaces; ++iFlx)
      {
        m_flxPntGradPtrs[iSide][iFlx] = &m_flxPntGrads[iSide][iFlx];
      }
    }
  }

  // set maximum number of flux points that will be passed at one time to the physical model
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
    setMaxNbrStatesData
    (
      prevMaxNbrStatesData > maxNbrFlxPntsFacesX2P2 ? prevMaxNbrStatesData : maxNbrFlxPntsFacesX2P2
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

void QuadFreeFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSide = 0; iSide < m_flxPntSol.size(); ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < m_flxPntSol[iSide].size(); ++iFlx)
    {
      deletePtr(m_flxPntSol      [iSide][iFlx]);
      deletePtr(m_flxPntExtraVars[iSide][iFlx]);
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
    for (CFuint iFlx = 0; iFlx < m_flxPntGrads[iSide].size(); ++iFlx)
    {
      for (CFuint iGrad = 0; iGrad < m_flxPntGrads[iSide][iFlx].size(); ++iGrad)
      {
        deletePtr(m_flxPntGrads[iSide][iFlx][iGrad]);
      }
      m_flxPntGrads[iSide][iFlx].resize(0);
    }
    m_flxPntGrads[iSide].resize(0);
  }
  m_flxPntGrads.resize(0);

  BaseFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
