#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/StdVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    StdVolTermComputer,SpectralFVMethodData,BaseVolTermComputer,SpectralFVModule >
  StdVolTermComputerProvider("StdVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

StdVolTermComputer::StdVolTermComputer(const std::string& name) :
  BaseVolTermComputer(name),
  m_localIntFaceCVConn(),
  m_intFaceQuadPntNorm(),
  m_intFaceQuadWheights(),
  m_normal(),
  m_faceFlux(),
  m_gradTerm()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

StdVolTermComputer::~StdVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void StdVolTermComputer::setVolumeTermData(CFuint iElemType)
{
  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // local internal face - CV connectivity
  m_localIntFaceCVConn = svLocalData[iElemType]->getLocalIntFaceCVConn();

  // local internal face normals in a SV
  m_intFaceQuadPntNorm = svLocalData[iElemType]->getLocalIntFaceNorm();

  // local internal face quadrature wheights
  m_intFaceQuadWheights = svLocalData[iElemType]->getIntFaceQuadWheights();

  // local internal quadrature quadrature point polynomial values
  m_flxPntsRecCoefs = svLocalData[iElemType]->getIntQuadPntPolyVals();

  // number of flux points in this element
  m_nbrFlxPnts = m_flxPntsRecCoefs->size();
}

//////////////////////////////////////////////////////////////////////////////

void StdVolTermComputer::computeConvVolTermFromFlxPntSol(RealVector& resUpdates)
{
  // set updates to zero
  resUpdates = 0.0;

  // compute fluxes through internal faces
  CFuint iState = 0;
  const CFuint nbrIntFaces = m_intFaceQuadWheights->size();
  for (CFuint iFace = 0; iFace < nbrIntFaces; ++iFace)
  {
    // compute flux integrated over the face
    m_faceFlux = 0.;

    // loop over quadrature points
    const CFuint nbrIntFaceQPnts = (*m_intFaceQuadWheights)[iFace].size();
    for (CFuint iQPnt = 0; iQPnt < nbrIntFaceQPnts; ++iQPnt, ++iState)
    {
      // compute local face normal
      m_normal = m_cellFaceNormTransfM * (*m_intFaceQuadPntNorm)[iFace][iQPnt];

      // dereference the state
      State& solQuadPnt = *m_solInFlxPnts[iState];

      // add contribution to total flux through face
      m_faceFlux += (*m_intFaceQuadWheights)[iFace][iQPnt]*m_updateVarSet->getFlux()(solQuadPnt,m_normal);
    }

    // add to the residual updates
    CFuint resIDL = m_nbrEqs*(*m_localIntFaceCVConn)[iFace][LEFT ];
    CFuint resIDR = m_nbrEqs*(*m_localIntFaceCVConn)[iFace][RIGHT];
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resIDL, ++resIDR)
    {
      resUpdates[resIDL] -= m_faceFlux[iVar];
      resUpdates[resIDR] += m_faceFlux[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdVolTermComputer::computeGradVolTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
{
  // set gradient volume terms to zero
  const CFuint nbrCVs = gradUpdates.size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradUpdates[iCV][iGrad] = 0.0;
    }
  }

  // compute gradient volume terms
  CFuint iState = 0;
  const CFuint nbrIntFaces = m_intFaceQuadWheights->size();
  for (CFuint iFace = 0; iFace < nbrIntFaces; ++iFace)
  {
    // left and right CVs
    const CFuint cvL = (*m_localIntFaceCVConn)[iFace][LEFT ];
    const CFuint cvR = (*m_localIntFaceCVConn)[iFace][RIGHT];

    // loop over quadrature points
    const CFuint nbrIntFaceQPnts = (*m_intFaceQuadWheights)[iFace].size();
    for (CFuint iQPnt = 0; iQPnt < nbrIntFaceQPnts; ++iQPnt, ++iState)
    {
      // compute local face normal
      m_normal = m_cellFaceNormTransfM * (*m_intFaceQuadPntNorm)[iFace][iQPnt];

      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
	// add contribution to gradient term
	m_gradTerm = (*m_intFaceQuadWheights)[iFace][iQPnt]*
	  m_gradVarsInFlxPnts(iGrad,iState)*m_normal;

        gradUpdates[cvL][iGrad] += m_gradTerm;
        gradUpdates[cvR][iGrad] -= m_gradTerm;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdVolTermComputer::computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates)
{
  // set updates to zero
  resUpdates = 0.0;

  // compute fluxes through internal faces
  CFuint iState = 0;
  const CFuint nbrIntFaces = m_intFaceQuadWheights->size();
  for (CFuint iFace = 0; iFace < nbrIntFaces; ++iFace)
  {
    // compute flux integrated over the face
    m_faceFlux = 0.;

    // loop over quadrature points
    const CFuint nbrIntFaceQPnts = (*m_intFaceQuadWheights)[iFace].size();
    for (CFuint iQPnt = 0; iQPnt < nbrIntFaceQPnts; ++iQPnt, ++iState)
    {
      // compute local face normal
      m_normal = m_cellFaceNormTransfM * (*m_intFaceQuadPntNorm)[iFace][iQPnt];

      // dereference the state
      RealVector& solQuadPnt = *m_solInFlxPnts[iState];

      // dereference the gradients
      vector< RealVector* >& gradQuadPnt = m_gradInFlxPnts[iState];

      // call the setComposition function of the DiffusiveVarSet
//      m_diffusiveVarSet->setComposition(...,false,0);// only needed for chemically reacting flows

      // add contribution to diffusive flux through face
      m_faceFlux += (*m_intFaceQuadWheights)[iFace][iQPnt]*
                                            m_diffusiveVarSet->getFlux(solQuadPnt,gradQuadPnt,m_normal,0);
    }

    // add to the residual updates
    CFuint resIDL = m_nbrEqs*(*m_localIntFaceCVConn)[iFace][LEFT ];
    CFuint resIDR = m_nbrEqs*(*m_localIntFaceCVConn)[iFace][RIGHT];
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resIDL, ++resIDR)
    {
      resUpdates[resIDL] += m_faceFlux[iVar];
      resUpdates[resIDR] -= m_faceFlux[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdVolTermComputer::setup()
{
  CFAUTOTRACE;

  BaseVolTermComputer::setup();

  // dimensionality and number of variables
  const CFuint dim     = PhysicalModelStack::getActive()->getDim ();

  // resize variables
  m_normal.resize(dim);
  m_faceFlux.resize(m_nbrEqs);

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  const CFuint nbrElemTypes = svLocalData.size();
  cf_assert(nbrElemTypes > 0);

  // get the maximum number of quadrature nodes
  CFuint maxNbrQNodes = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    m_flxPntsRecCoefs = svLocalData[iElemType]->getIntQuadPntPolyVals();
    const CFuint nbrQNodes = m_flxPntsRecCoefs->size();
    maxNbrQNodes = maxNbrQNodes > nbrQNodes ? maxNbrQNodes : nbrQNodes;
  }

  // create states and extra variables for quadrature points
  m_solInFlxPnts      .resize(maxNbrQNodes);
  m_extraVarsInFlxPnts.resize(maxNbrQNodes);
  for (CFuint iState = 0; iState < maxNbrQNodes; ++iState)
  {
    m_solInFlxPnts      [iState] = new State();
    m_extraVarsInFlxPnts[iState] = new RealVector(m_nbrExtraVars);
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(maxNbrQNodes);

  // setup variables for gradient computation
  if (getMethodData().hasDiffTerm())
  {
    // create gradient vars for quadrature points
    m_gradVarsInFlxPnts.resize(m_nbrEqs, maxNbrQNodes);

    // set RealVector pointers to states
    m_solRVInFlxPnts.resize(maxNbrQNodes);
    for (CFuint iQPnt = 0; iQPnt < maxNbrQNodes; ++iQPnt)
    {
      m_solRVInFlxPnts[iQPnt] = m_solInFlxPnts[iQPnt];
    }

    // create gradients for quadrature points
    m_gradInFlxPnts.resize(maxNbrQNodes);
    for (CFuint iQPnt = 0; iQPnt < maxNbrQNodes; ++iQPnt)
    {
      m_gradInFlxPnts[iQPnt].resize(m_nbrEqs);
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        m_gradInFlxPnts[iQPnt][iGrad] = new RealVector(dim);
      }
    }

    // resize m_gradTerm
    m_gradTerm.resize(dim);
  }

  // set maximum number of quadrature points that will be passed at one time to the physical model
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
    setMaxNbrStatesData(prevMaxNbrStatesData > maxNbrQNodes ? prevMaxNbrStatesData : maxNbrQNodes);
}

//////////////////////////////////////////////////////////////////////////////

void StdVolTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iState = 0; iState < m_solInFlxPnts.size(); ++iState)
  {
    deletePtr(m_solInFlxPnts      [iState]);
    deletePtr(m_extraVarsInFlxPnts[iState]);
  }
  m_solInFlxPnts      .resize(0);
  m_extraVarsInFlxPnts.resize(0);

  for (CFuint iQPnt = 0; iQPnt < m_gradInFlxPnts.size(); ++iQPnt)
  {
    for (CFuint iGrad = 0; iGrad < m_gradInFlxPnts[iQPnt].size(); ++iGrad)
    {
      deletePtr(m_gradInFlxPnts[iQPnt][iGrad]);
    }
    m_gradInFlxPnts[iQPnt].resize(0);
  }
  m_gradInFlxPnts.resize(0);

  BaseVolTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

