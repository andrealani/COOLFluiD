#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/CompactVolTermComputer.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
  CompactVolTermComputer,SpectralFDMethodData,BaseVolTermComputer,SpectralFDModule >
  CompactVolTermComputerProvider("CompactVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

CompactVolTermComputer::CompactVolTermComputer(const std::string& name) :
  BaseVolTermComputer(name),
  m_backupPhysVarGrad()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CompactVolTermComputer::~CompactVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CompactVolTermComputer::backupAndReconstructGradients(const CFuint iVar,
                                                           const vector< vector< RealVector >* >& cellGradients)
{
  // backup
  backupPhysVarGrad(iVar);

  // reconstruct
  m_statesReconstr->reconstructPhysVarGrad(iVar,cellGradients,m_gradInFlxPnts,
                                           *m_flxPntsRecCoefs,*m_intFlxPntIdxs,
                                           *m_flxPntMatrixIdxForReconstruction,
                                           *m_solPntIdxsForReconstruction);
}

//////////////////////////////////////////////////////////////////////////////

void CompactVolTermComputer::backupPhysVarGrad(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts <= m_backupPhysVarGrad.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    m_backupPhysVarGrad[iFlx] = *m_gradInFlxPnts[iFlx][iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactVolTermComputer::restorePhysVarGrad(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts <= m_backupPhysVarGrad.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    *m_gradInFlxPnts[iFlx][iVar] = m_backupPhysVarGrad[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactVolTermComputer::computePhysVarGradientVolumeTerm(const CFuint iVar, vector< RealVector >& physVarGradUpdates)
{
  computePhysVarGradVolTermFromFlxPntSol(iVar,physVarGradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void CompactVolTermComputer::computePhysVarGradVolTermFromFlxPntSol(const CFuint iVar, vector< RealVector >& physVarGradUpdates)
{
  // set gradient volume terms to zero
  const CFuint nbrSolPnts = physVarGradUpdates.size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    physVarGradUpdates[iSol] = 0.0;
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

    // compute gradient term
    m_gradTerm = (*m_solRVInFlxPnts[iFlx])[iVar]*m_cellFluxProjVects[iFlx];

    // add contribution of this flux point to the solution points
    for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
    {
      // solution point index
      const CFuint solIdx = (*m_solPntIdxsForDerivation)[flxIdx][iSol];
      physVarGradUpdates[solIdx] += (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*m_gradTerm;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactVolTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  BaseVolTermComputer::setup();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  const CFuint nbrElemTypes = sdLocalData.size();
  cf_assert(nbrElemTypes > 0);

  // get the maximum number of flux points
  CFuint maxNbrFlxPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrFlxPnts = sdLocalData[iElemType]->getNbrOfIntFlxPnts();
    maxNbrFlxPnts = maxNbrFlxPnts > nbrFlxPnts ? maxNbrFlxPnts : nbrFlxPnts;
  }

  // resize m_backupPhysVarGrad
  m_backupPhysVarGrad.resize(maxNbrFlxPnts,RealVector(m_dim));
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
