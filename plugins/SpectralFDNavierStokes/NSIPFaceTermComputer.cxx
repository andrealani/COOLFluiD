#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NSIPFaceTermComputer.hh"

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
    NSIPFaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDNavierStokesModule >
NSIPFaceTermComputerProvider("NSIPFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSIPFaceTermComputer::NSIPFaceTermComputer(const std::string& name) :
  IPFaceTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSIPFaceTermComputer::~NSIPFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSIPFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(vector< RealVector>& resUpdates,
                                                                              vector< CFreal >& updateCoefContr)
{
  // compute the face term
  computeDiffFaceTerm(resUpdates);

  // compute the update coefficient contributions for the neighbouring cells
  cf_assert(updateCoefContr.size() == 2);
  const CFreal dynVisc = m_navierStokesVarSet->getCurrDynViscosity();
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    updateCoefContr[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         (*m_faceIntegrationCoefs)[iFlx]*m_cflConvDiffRatio;
      const CFreal rho = m_navierStokesVarSet->getDensity(*m_flxPntSol[iSide][iFlx]);
      updateCoefContr[iSide] += dynVisc*jacobXJacobXIntCoef/rho/m_cellVolumes[iSide];
      /// @todo replace cell volume with Jacobian determinant at the flux point
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSIPFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  IPFaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void NSIPFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of the parent class
  IPFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
