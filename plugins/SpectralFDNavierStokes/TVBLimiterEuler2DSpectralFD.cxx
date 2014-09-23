#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/TVBLimiterEuler2DSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TVBLimiterEuler2DSpectralFD, SpectralFDMethodData, SpectralFDNavierStokesModule>
    TVBLimiterEuler2DSpectralFDProvider("TVBLimiterEuler2D");

//////////////////////////////////////////////////////////////////////////////

TVBLimiterEuler2DSpectralFD::TVBLimiterEuler2DSpectralFD(const std::string& name) :
  TVBLimiterEulerSpectralFD(name),
  m_eulerVarSet(CFNULL),
  m_gammaMinusOne()
{
}

//////////////////////////////////////////////////////////////////////////////

TVBLimiterEuler2DSpectralFD::~TVBLimiterEuler2DSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2DSpectralFD::configure ( Config::ConfigArgs& args )
{
  TVBLimiterSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2DSpectralFD::computePressFromConsVar(const RealVector& consVar, CFreal& press)
{
  const CFreal& rho  = consVar[0];
  const CFreal& rhoU = consVar[1];
  const CFreal& rhoV = consVar[2];
  const CFreal& rhoE = consVar[3];

  press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2DSpectralFD::computeRhoEFromPressAndOtherConsVar(RealVector& consVar, const CFreal& press)
{
  const CFreal& rho  = consVar[0];
  const CFreal& rhoU = consVar[1];
  const CFreal& rhoV = consVar[2];
  CFreal& rhoE = consVar[3];

  rhoE = press/m_gammaMinusOne + 0.5*(rhoU*rhoU+rhoV*rhoV)/rho;
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2DSpectralFD::setup()
{
  CFAUTOTRACE;

  TVBLimiterEulerSpectralFD::setup();

  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in TVBLimiterEuler2DSpectralFD!");
  }
  cf_assert(m_nbrEqs == 4);

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2DSpectralFD::unsetup()
{
  CFAUTOTRACE;

  TVBLimiterEulerSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
