#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/NavierStokesVolTermComputer.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/LESVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LESVolTermComputer,SpectralFDMethodData,BaseVolTermComputer,SpectralFDLESModule >
LESVolTermComputerProvider("LESVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

LESVolTermComputer::LESVolTermComputer(const std::string& name) :
  NavierStokesVolTermComputer(name),
  m_lesVarSet(CFNULL),
  m_nbrOfIntFlxPnts(),
  m_filterWidthVolumes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LESVolTermComputer::~LESVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LESVolTermComputer::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NavierStokesVolTermComputer::setup();

  // Dynamic_cast the diffusive varset to an LES varset
  m_lesVarSet = m_diffusiveVarSet.d_castTo< LES::LESVarSet >();

}

//////////////////////////////////////////////////////////////////////////////

void LESVolTermComputer::setVolumeTermData(CFuint iElemType)
{
  CFAUTOTRACE;

  // Call setVolumeTermData() of the parent class
  NavierStokesVolTermComputer::setVolumeTermData(iElemType);

  // Get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // Number of internal flux points in this element
  m_nbrOfIntFlxPnts = static_cast<CFreal>(sdLocalData[iElemType]->getNbrOfIntFlxPnts());

  // Resize m_filterWidthVolumes
  m_filterWidthVolumes.resize(m_nbrOfIntFlxPnts);

}

//////////////////////////////////////////////////////////////////////////////

void LESVolTermComputer::computeCellData()
{
  CFAUTOTRACE;

  // Call computeCellData() of the parent class
  NavierStokesVolTermComputer::computeCellData();

  // Compute Jacobian determinants in the internal flux points
  std::valarray<CFreal> jacobDets(m_nbrFlxPnts);
  jacobDets = m_cell->computeGeometricShapeFunctionJacobianDeterminant(m_intFlxPntMappedCoord);

  // set averaged jacobian determinants in the LESVarSet
  for (CFuint iFlx = 0; iFlx < m_nbrOfIntFlxPnts; ++iFlx)
  {
    m_filterWidthVolumes[iFlx] = jacobDets[iFlx];
  }

}

//////////////////////////////////////////////////////////////////////////////

void LESVolTermComputer::computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates)
{
  CFAUTOTRACE;

  // Set updates to zero
  resUpdates = 0.0;

  // Number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // Loop over flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    // Flux point index
    const CFuint flxIdx = (*m_intFlxPntIdxs)[iFlx];

    // Flux point index in the matrix m_solPntsDerivCoefs
    const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

    // Dereference the state in the flux point
    State& solFlxPnt = *m_solInFlxPnts[iFlx];

    // Dereference the gradients
    vector< RealVector* >& gradFlxPnt = m_gradVarGradsInFlxPnts[iFlx];

    // Filter Width Volume to use in the LES calculation
    m_lesVarSet->setVolume(m_filterWidthVolumes[iFlx]);

    // Evaluate flux projected on the projection vector
    const RealVector fluxXProjVect
        = m_diffusiveVarSet->getFlux(solFlxPnt,gradFlxPnt,m_cellFluxProjVects[iFlx],0.);

    // Add contribution of this flux point to the solution points
    for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
    {
      // First residual in solution point index
      CFuint resIdx = m_nbrEqs*(*m_solPntIdxsForDerivation)[flxIdx][iSol];
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resIdx)
      {
        resUpdates[resIdx] += (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*fluxXProjVect[iVar];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LESVolTermComputer::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NavierStokesVolTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
