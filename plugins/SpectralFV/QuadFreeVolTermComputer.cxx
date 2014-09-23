#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/QuadFreeVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    QuadFreeVolTermComputer,SpectralFVMethodData,BaseVolTermComputer,SpectralFVModule >
  QuadFreeVolTermComputerProvider("QuadFreeVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

QuadFreeVolTermComputer::QuadFreeVolTermComputer(const std::string& name) :
  BaseVolTermComputer(name),
  m_volTermTensor(),
  m_physFlux(),
  m_faceFlux(),
  m_dim()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

QuadFreeVolTermComputer::~QuadFreeVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeVolTermComputer::setVolumeTermData(CFuint iElemType)
{
  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // coefficients for solution reconstruction in the flux points
  m_flxPntsRecCoefs = svLocalData[iElemType]->getSolInFluxCoef();

    // tensor for the evaluation of the volume terms
  m_volTermTensor = svLocalData[iElemType]->getVolTermTensor();

  // number of flux points in this element
  m_nbrFlxPnts = m_flxPntsRecCoefs->size();
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeVolTermComputer::computeConvVolTermFromFlxPntSol(RealVector& resUpdates)
{
  // set updates to zero
  resUpdates = 0.0;

  // compute volume term
  const CFuint nbrCVs = m_volTermTensor->size();
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    // dereference the state
    State& solInFlxPnt = *m_solInFlxPnts[iFlx];

    // compute physical flux from solution in flux point
    m_physFlux = m_updateVarSet->getFlux()(solInFlxPnt);

    // compute mapped flux coefficients
    for (CFuint iMappedCoor = 0; iMappedCoor < m_dim; ++iMappedCoor)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        CFreal mappedFluxCoef = 0.0;
        for (CFuint iAbsCoor = 0; iAbsCoor < m_dim; ++iAbsCoor)
        {
          mappedFluxCoef += m_cellFaceNormTransfM(iAbsCoor,iMappedCoor)*m_physFlux(iVar,iAbsCoor);
        }

        for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
        {
          const CFuint resID = m_nbrEqs*iCV + iVar;
          resUpdates[resID] += (*m_volTermTensor)[iCV][iFlx][iMappedCoor]*mappedFluxCoef;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeVolTermComputer::computeGradVolTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
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

  // add volume terms to gradients
  cf_assert(m_gradVarsInFlxPnts.nbRows() > 0);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
    {
      for (CFuint iCoor = 0; iCoor < m_dim; ++iCoor)
      {
        CFreal coef = 0.0;
        for (CFuint iMappedCoor = 0; iMappedCoor < m_dim; ++iMappedCoor)
        {
          coef += (*m_volTermTensor)[iCV][iFlx][iMappedCoor]*m_cellFaceNormTransfM(iCoor,iMappedCoor);
        }

        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
        {
          // add contribution to gradient
          gradUpdates[iCV][iGrad][iCoor] -= coef*m_gradVarsInFlxPnts(iGrad,iFlx);
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeVolTermComputer::computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates)
{
  // get the dimensionality
  const CFuint m_dim = PhysicalModelStack::getActive()->getDim();

  // set updates to zero
  resUpdates = 0.0;

  // compute volume term
  const CFuint nbrCVs = m_volTermTensor->size();
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    // dereference the state
    RealVector& solInFlxPnt = *m_solInFlxPnts[iFlx];

    // dereference the gradients
    vector< RealVector* >& gradInFlxPnt = m_gradInFlxPnts[iFlx];

    // call the setComposition function of the DiffusiveVarSet
//    m_diffusiveVarSet->setComposition(...,false,0);// only needed for chemically reacting flows

    // compute physical flux from solution in flux point
    m_physFlux = m_diffusiveVarSet->getFlux(solInFlxPnt,gradInFlxPnt,0);

    // compute mapped flux coefficients
    for (CFuint iMappedCoor = 0; iMappedCoor < m_dim; ++iMappedCoor)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        CFreal mappedFluxCoef = 0.0;
        for (CFuint iAbsCoor = 0; iAbsCoor < m_dim; ++iAbsCoor)
        {
          mappedFluxCoef += m_cellFaceNormTransfM(iAbsCoor,iMappedCoor)*m_physFlux(iVar,iAbsCoor);
        }

        for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
        {
          const CFuint resID = m_nbrEqs*iCV + iVar;
          resUpdates[resID] -= (*m_volTermTensor)[iCV][iFlx][iMappedCoor]*mappedFluxCoef;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeVolTermComputer::setup()
{
  CFAUTOTRACE;

  BaseVolTermComputer::setup();

  // dimensionality and number of variables
  m_dim = PhysicalModelStack::getActive()->getDim ();

  // resize variables
  m_physFlux.resize(m_nbrEqs,m_dim);
  m_faceFlux.resize(m_nbrEqs);

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  const CFuint nbrElemTypes = svLocalData.size();
  cf_assert(nbrElemTypes > 0);

  // maximum number of flux points
  CFuint maxNbrFlxPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get number of flux points
    const CFuint nbrFlxPnts = svLocalData[iElemType]->getNbrOfFlxPnts();
    maxNbrFlxPnts = maxNbrFlxPnts > nbrFlxPnts ? maxNbrFlxPnts : nbrFlxPnts;
  }

  // create states for flux points
  m_solInFlxPnts      .resize(maxNbrFlxPnts);
  m_extraVarsInFlxPnts.resize(maxNbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPnts; ++iFlx)
  {
    m_solInFlxPnts      [iFlx] = new State();
    m_extraVarsInFlxPnts[iFlx] = new RealVector(m_nbrExtraVars);
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(maxNbrFlxPnts);

  // setup variables for gradient computation
  if (getMethodData().hasDiffTerm())
  {
    // create gradient vars for flux points
    m_gradVarsInFlxPnts.resize(m_nbrEqs,maxNbrFlxPnts);

    // set RealVector pointers to states
    m_solRVInFlxPnts.resize(maxNbrFlxPnts);
    for (CFuint iFlx = 0; iFlx < maxNbrFlxPnts; ++iFlx)
    {
      m_solRVInFlxPnts[iFlx] = m_solInFlxPnts[iFlx];
    }

    // create gradients for quadrature points
    m_gradInFlxPnts.resize(maxNbrFlxPnts);
    for (CFuint iFlx = 0; iFlx < maxNbrFlxPnts; ++iFlx)
    {
      m_gradInFlxPnts[iFlx].resize(m_nbrEqs);
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        m_gradInFlxPnts[iFlx][iGrad] = new RealVector(m_dim);
      }
    }
  }

  // set maximum number of quadrature points that will be passed at one time to the physical model
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
    setMaxNbrStatesData(prevMaxNbrStatesData > maxNbrFlxPnts ? prevMaxNbrStatesData : maxNbrFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void QuadFreeVolTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iState = 0; iState < m_solInFlxPnts.size(); ++iState)
  {
    deletePtr(m_solInFlxPnts      [iState]);
    deletePtr(m_extraVarsInFlxPnts[iState]);
  }
  m_solInFlxPnts      .resize(0);
  m_extraVarsInFlxPnts.resize(0);

  for (CFuint iFlx = 0; iFlx < m_gradInFlxPnts.size(); ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_gradInFlxPnts[iFlx].size(); ++iGrad)
    {
      deletePtr(m_gradInFlxPnts[iFlx][iGrad]);
    }
    m_gradInFlxPnts[iFlx].resize(0);
  }
  m_gradInFlxPnts.resize(0);

  BaseVolTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

