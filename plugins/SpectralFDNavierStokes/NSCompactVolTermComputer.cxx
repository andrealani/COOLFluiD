#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NSCompactVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    NSCompactVolTermComputer,SpectralFDMethodData,BaseVolTermComputer,SpectralFDNavierStokesModule >
NSCompactVolTermComputerProvider("NSCompactVolTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSCompactVolTermComputer::NSCompactVolTermComputer(const std::string& name) :
  CompactVolTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSCompactVolTermComputer::~NSCompactVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSCompactVolTermComputer::computeCellDiffVolumeTerm(RealVector& resUpdates)
{
  // transform the cell state gradients to gradient variable gradients
  m_navierStokesVarSet->setGradientVarGradients(m_solRVInFlxPnts,
                                                m_gradInFlxPnts,
                                                m_gradVarGradsInFlxPnts,
                                                m_nbrFlxPnts);

  // compute the actual volume term
  computeDiffVolTermFromFlxPntSolAndGrad(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void NSCompactVolTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  CompactVolTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

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

  // create gradients for flux points
  for (CFuint iFlx = 0; iFlx < maxNbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_gradVarGradsInFlxPnts[iFlx][iGrad] = new RealVector(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSCompactVolTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iFlx = 0; iFlx < m_gradVarGradsInFlxPnts.size(); ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_gradVarGradsInFlxPnts[iFlx].size(); ++iGrad)
    {
      deletePtr(m_gradVarGradsInFlxPnts[iFlx][iGrad]);
    }
    m_gradVarGradsInFlxPnts[iFlx].resize(0);
  }
  m_gradVarGradsInFlxPnts.resize(0);

  // call unsetup of the parent class
  CompactVolTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
