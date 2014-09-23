#include "Common/BadValueException.hh"
#include "Framework/ConsistencyException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BaseTerm.hh"
#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BaseVolTermComputer,SpectralFDMethodData,BaseVolTermComputer,SpectralFDModule >
    BaseVolTermComputerProvider("BaseVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

BaseVolTermComputer::BaseVolTermComputer(const std::string& name) :
  SpectralFDMethodStrategy(name),
  socket_extraVars("meanflow",false),/// @todo get the name of this socket in another way, not hardcoded, to make it more general
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_statesReconstr(CFNULL),
  m_flxPntsRecCoefs(CFNULL),
  m_solPntsDerivCoefs(CFNULL),
  m_intFlxPntIdxs(CFNULL),
  m_flxPntDerivDir(CFNULL),
  m_intFlxPntDerivDir(),
  m_flxPntsLocalCoords(),
  m_intFlxPntMappedCoord(),
  m_flxPntMatrixIdxForReconstruction(CFNULL),
  m_solPntIdxsForReconstruction(CFNULL),
  m_flxPntMatrixIdxForDerivation(CFNULL),
  m_solPntIdxsForDerivation(CFNULL),
  m_cellFluxProjVects(),
  m_cellExtraVars(),
  m_solInFlxPnts(),
  m_solRVInFlxPnts(),
  m_extraVarsInFlxPnts(),
  m_gradVarsInFlxPnts(),
  m_gradInFlxPnts(),
  m_gradVarGradsInFlxPnts(),
  m_backupPhysVar(),
  m_nbrFlxPnts(),
  m_nbrEqs(),
  m_dim(),
  m_gradTerm(),
  m_nbrExtraVars()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BaseVolTermComputer::~BaseVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodStrategy::configure(args);

  // initialize socket_extraVars
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::setVolumeTermData(CFuint iElemType)
{
  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get derivation coefficients for the solution points
  m_solPntsDerivCoefs = sdLocalData[iElemType]->getDerivCoefsSolPnts1D();

  // get indexes of internal flux points
  m_intFlxPntIdxs = sdLocalData[iElemType]->getIntFlxPntIdxs();

  // get derivation direction in the flux points
  m_flxPntDerivDir = sdLocalData[iElemType]->getFlxPntDerivDir();

  // get flux points mapped coordinates
  m_flxPntsLocalCoords = sdLocalData[iElemType]->getFlxPntsLocalCoords();

  // get flux point index (in the matrix flxPntRecCoefs) for reconstruction
  m_flxPntMatrixIdxForReconstruction = sdLocalData[iElemType]->getFlxPntMatrixIdxForReconstruction();

  // get flux point index (in the matrix m_solPntsDerivCoefs) for derivation
  m_flxPntMatrixIdxForDerivation = sdLocalData[iElemType]->getFlxPntMatrixIdxForDerivation();

  // get solution point index (in the cell) for derivation
  m_solPntIdxsForDerivation = sdLocalData[iElemType]->getSolPntIdxsForDerivation();

  // number of internal flux points in this element
  m_nbrFlxPnts = sdLocalData[iElemType]->getNbrOfIntFlxPnts();

  // set derivation direction and mapped coordinates in the internal flux points
  m_intFlxPntDerivDir.resize(0);
  m_intFlxPntMappedCoord.resize(0);
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_intFlxPntIdxs)[iFlx];
    m_intFlxPntDerivDir.push_back((*m_flxPntDerivDir)[flxIdx]);
    m_intFlxPntMappedCoord.push_back((*m_flxPntsLocalCoords)[flxIdx]);
  }

  // set interpolation type dependent data
  const std::string interpolationType = getMethodData().getInterpolationType();
  if (interpolationType == "standard")
  {
    // get reconstruction coefficients for the flux points
    m_flxPntsRecCoefs = sdLocalData[iElemType]->getRecCoefsFlxPnts1D();

    // get solution point index (in the cell) for reconstruction
    m_solPntIdxsForReconstruction = sdLocalData[iElemType]->getSolPntIdxsForReconstruction();
  }
  else if (interpolationType == "optimized")
  {
    // get reconstruction coefficients for the flux points
    m_flxPntsRecCoefs = sdLocalData[iElemType]->getRecCoefsFlxPnts1DOptim();

    // get solution point index (in the cell) for reconstruction
    m_solPntIdxsForReconstruction = sdLocalData[iElemType]->getSolPntIdxsForRecOptim();
  }
  else
  {
    throw BadValueException (FromHere(),"BaseVolTermComputer::setVolumeTermData --> Interpolation type should be standard or optimized");
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeCellData()
{
  m_cellFluxProjVects = m_cell->computeMappedCoordPlaneNormalAtMappedCoords(m_intFlxPntDerivDir,
                                                                            m_intFlxPntMappedCoord);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::reconstructStates(const vector< State* >& cellStates, bool onlyExtraVars)
{
  if (!onlyExtraVars)
  {
    m_statesReconstr->reconstructStates(cellStates,m_solInFlxPnts,
                                        *m_flxPntsRecCoefs,*m_intFlxPntIdxs,
                                        *m_flxPntMatrixIdxForReconstruction,
                                        *m_solPntIdxsForReconstruction);
  }

  // if needed, reconstruct the extra variables
  if (m_nbrExtraVars > 0)
  {
    cf_assert(socket_extraVars.isConnected());
    DataHandle<RealVector> extraVars = socket_extraVars.getDataHandle();

    // get extra vars at the solution points in this cell
    const CFuint nbrSolPnts = cellStates.size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFuint stateID = cellStates[iSol]->getLocalID();
      m_cellExtraVars[iSol] = &extraVars[stateID];
    }

    // reconstruct extra variables
    m_statesReconstr->reconstructExtraVars(m_cellExtraVars,m_extraVarsInFlxPnts,
                                           *m_flxPntsRecCoefs,*m_intFlxPntIdxs,
                                           *m_flxPntMatrixIdxForReconstruction,
                                           *m_solPntIdxsForReconstruction);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::reconstructGradients(const vector< vector< RealVector >* >& cellGradients)
{
  m_statesReconstr->reconstructGradients(cellGradients,m_gradInFlxPnts,
                                         *m_flxPntsRecCoefs,*m_intFlxPntIdxs,
                                         *m_flxPntMatrixIdxForReconstruction,
                                         *m_solPntIdxsForReconstruction);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::backupAndReconstructPhysVar(const CFuint iVar, const vector< State* >& cellStates)
{
  // backup
  backupPhysVar(iVar);

  // reconstruct
  m_statesReconstr->reconstructPhysVar(iVar,cellStates,m_solInFlxPnts,
                                       *m_flxPntsRecCoefs,*m_intFlxPntIdxs,
                                       *m_flxPntMatrixIdxForReconstruction,
                                       *m_solPntIdxsForReconstruction);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::backupPhysVar(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    m_backupPhysVar[iFlx] = (*m_solInFlxPnts[iFlx])[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::restorePhysVar(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    (*m_solInFlxPnts[iFlx])[iVar] = m_backupPhysVar[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeCellConvVolumeTerm(RealVector& resUpdates)
{
  // compute the actual volume term
  computeConvVolTermFromFlxPntSol(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeGradientVolumeTerm(vector< vector< RealVector > >& gradUpdates)
{
  // set the gradient variables in the flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_gradVarsInFlxPnts(iGrad,iFlx) = (*m_solRVInFlxPnts[iFlx])[iGrad];
    }
  }

  // compute the actual volume term contribution to the gradient
  computeGradVolTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeGradientExtraVarsVolumeTerm(vector< vector< RealVector > >& gradUpdates)
{
  // set the extra variables in the flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_gradVarsInFlxPnts(iGrad,iFlx) = (*m_extraVarsInFlxPnts[iFlx])[iGrad];
    }
  }

  // compute the actual volume term contribution to the extra variable gradients
  computeGradVolTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeCellDiffVolumeTerm(RealVector& resUpdates)
{
  // compute the actual volume term
  computeDiffVolTermFromFlxPntSolAndGrad(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeConvVolTermFromFlxPntSol(RealVector& resUpdates)
{
  // set updates to zero
  resUpdates = 0.0;

  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // loop over flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    // flux point index
    const CFuint flxIdx = (*m_intFlxPntIdxs)[iFlx];

    // flux point index in the matrix m_solPntsDerivCoefs
    const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

    // dereference the state in the flux point
    State& solFlxPnt = *m_solInFlxPnts[iFlx];

    // transform update states to physical data to calculate eigenvalues
    m_updateVarSet->computePhysicalData(solFlxPnt, m_pData);

    // evaluate flux projected on the projection vector
    const RealVector fluxXProjVect
        = m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVects[iFlx]);

    // add contribution of this flux point to the solution points
    for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
    {
      // first residual in solution point index
      CFuint resIdx = m_nbrEqs*(*m_solPntIdxsForDerivation)[flxIdx][iSol];
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resIdx)
      {
        resUpdates[resIdx] -= (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*fluxXProjVect[iVar];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeGradVolTermFromFlxPntSol(vector< vector< RealVector > >& gradUpdates)
{
  // set gradient volume terms to zero
  const CFuint nbrSolPnts = gradUpdates.size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradUpdates[iSol][iGrad] = 0.0;
    }
  }

  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // loop over flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    // flux point index
    const CFuint flxIdx = (*m_intFlxPntIdxs)[iFlx];

    // flux point index in the matrix m_solPntsDerivCoefs
    const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute gradient term
      m_gradTerm = m_gradVarsInFlxPnts(iGrad,iFlx)*m_cellFluxProjVects[iFlx];

      // add contribution of this flux point to the solution points
      for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
      {
        // solution point index
        const CFuint solIdx = (*m_solPntIdxsForDerivation)[flxIdx][iSol];
        gradUpdates[solIdx][iGrad] += (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*m_gradTerm;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates)
{
  // set updates to zero
  resUpdates = 0.0;

  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // loop over flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    // flux point index
    const CFuint flxIdx = (*m_intFlxPntIdxs)[iFlx];

    // flux point index in the matrix m_solPntsDerivCoefs
    const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

    // dereference the state in the flux point
    State& solFlxPnt = *m_solInFlxPnts[iFlx];

      // dereference the gradients
    vector< RealVector* >& gradFlxPnt = m_gradVarGradsInFlxPnts[iFlx];

    // evaluate flux projected on the projection vector
    const RealVector fluxXProjVect
        = m_diffusiveVarSet->getFlux(solFlxPnt,gradFlxPnt,m_cellFluxProjVects[iFlx],0.);

    // add contribution of this flux point to the solution points
    for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
    {
      // first residual in solution point index
      CFuint resIdx = m_nbrEqs*(*m_solPntIdxsForDerivation)[flxIdx][iSol];
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resIdx)
      {
        resUpdates[resIdx] += (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*fluxXProjVect[iVar];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::setup()
{
  CFAUTOTRACE;

  // dimensionality and number of variables
  m_dim    = PhysicalModelStack::getActive()->getDim ();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();

  // get number of extra variables from the update variable set
  m_nbrExtraVars = m_updateVarSet->getExtraPhysicalVarsSize();
  if (m_nbrExtraVars > 0 && !socket_extraVars.isConnected())
  {
    throw ConsistencyException(FromHere(), "Extra variables socket not connected but required...");
  }

  // get the states reconstructor
  m_statesReconstr = getMethodData().getStatesReconstructor();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  const CFuint nbrElemTypes = sdLocalData.size();
  cf_assert(nbrElemTypes > 0);

  // get the maximum number of flux points
  CFuint maxNbrSolPnts = 0;
  CFuint maxNbrFlxPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrSolPnts = sdLocalData[iElemType]->getNbrOfSolPnts();
    maxNbrSolPnts = maxNbrSolPnts > nbrSolPnts ? maxNbrSolPnts : nbrSolPnts;

    const CFuint nbrFlxPnts = sdLocalData[iElemType]->getNbrOfIntFlxPnts();
    maxNbrFlxPnts = maxNbrFlxPnts > nbrFlxPnts ? maxNbrFlxPnts : nbrFlxPnts;
  }

  // resize m_cellExtraVars
  m_cellExtraVars.resize(maxNbrSolPnts);

  // create states and extra variables for flux points
  for (CFuint iState = 0; iState < maxNbrFlxPnts; ++iState)
  {
    m_solInFlxPnts      .push_back(new State()                   );
    m_extraVarsInFlxPnts.push_back(new RealVector(m_nbrExtraVars));
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(maxNbrFlxPnts);

  // setup variables for gradient computation
  if (getMethodData().hasDiffTerm())
  {
    // get the diffusive varset
    m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  }

  // resize m_gradVarsInFlxPnts
  m_gradVarsInFlxPnts.resize(m_nbrEqs, maxNbrFlxPnts);

  // set RealVector pointers to states
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPnts; ++iFlx)
  {
    m_solRVInFlxPnts.push_back(m_solInFlxPnts[iFlx]);
  }

  // create gradients for flux points
  m_gradInFlxPnts        .resize(maxNbrFlxPnts);
  m_gradVarGradsInFlxPnts.resize(maxNbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_gradInFlxPnts        [iFlx].push_back(new RealVector(m_dim));
      m_gradVarGradsInFlxPnts[iFlx].push_back(m_gradInFlxPnts[iFlx][iGrad]);
    }
  }

  // resize m_gradTerm
  m_gradTerm.resize(m_dim);

  // set maximum number of states in flux points that will be passed at one time to the physical model
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
      setMaxNbrStatesData(prevMaxNbrStatesData > maxNbrFlxPnts ? prevMaxNbrStatesData : maxNbrFlxPnts);
  
  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::unsetup()
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
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > BaseVolTermComputer::needsSockets()
{
  vector< SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_extraVars);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
