#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NSBR2FaceTermComputer.hh"

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
    NSBR2FaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDNavierStokesModule >
NSBR2FaceTermComputerProvider("NSBR2FaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

NSBR2FaceTermComputer::NSBR2FaceTermComputer(const std::string& name) :
  BR2FaceTermComputer(name),
  m_navierStokesVarSet(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSBR2FaceTermComputer::~NSBR2FaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NSBR2FaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(vector< RealVector>& resUpdates,
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

void NSBR2FaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of the parent class
  BR2FaceTermComputer::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void NSBR2FaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of the parent class
  BR2FaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
