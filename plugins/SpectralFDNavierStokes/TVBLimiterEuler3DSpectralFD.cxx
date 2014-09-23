#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "NavierStokes/Euler3DVarSet.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/TVBLimiterEuler3DSpectralFD.hh"

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

MethodCommandProvider<TVBLimiterEuler3DSpectralFD, SpectralFDMethodData, SpectralFDNavierStokesModule>
    TVBLimiterEuler3DSpectralFDProvider("TVBLimiterEuler3D");

//////////////////////////////////////////////////////////////////////////////

TVBLimiterEuler3DSpectralFD::TVBLimiterEuler3DSpectralFD(const std::string& name) :
  TVBLimiterEulerSpectralFD(name),
  m_eulerVarSet(CFNULL),
  m_gammaMinusOne()
{
}

//////////////////////////////////////////////////////////////////////////////

TVBLimiterEuler3DSpectralFD::~TVBLimiterEuler3DSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler3DSpectralFD::configure ( Config::ConfigArgs& args )
{
  TVBLimiterSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler3DSpectralFD::computePressFromConsVar(const RealVector& consVar, CFreal& press)
{
  const CFreal& rho  = consVar[0];
  const CFreal& rhoU = consVar[1];
  const CFreal& rhoV = consVar[2];
  const CFreal& rhoW = consVar[3];
  const CFreal& rhoE = consVar[4];

  press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler3DSpectralFD::computeRhoEFromPressAndOtherConsVar(RealVector& consVar, const CFreal& press)
{
  const CFreal& rho  = consVar[0];
  const CFreal& rhoU = consVar[1];
  const CFreal& rhoV = consVar[2];
  const CFreal& rhoW = consVar[3];
  CFreal& rhoE = consVar[4];

  rhoE = press/m_gammaMinusOne + 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho;
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler3DSpectralFD::setup()
{
  CFAUTOTRACE;

  TVBLimiterEulerSpectralFD::setup();

  // get Euler 3D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in TVBLimiterEuler3DSpectralFD!");
  }
  cf_assert(m_nbrEqs == 4);

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler3DSpectralFD::unsetup()
{
  CFAUTOTRACE;

  TVBLimiterEulerSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
