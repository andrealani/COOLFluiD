#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/TVBLimiterEuler2D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TVBLimiterEuler2D, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    TVBLimiterEuler2DFRProvider("TVBLimiterEuler2D");

//////////////////////////////////////////////////////////////////////////////

TVBLimiterEuler2D::TVBLimiterEuler2D(const std::string& name) :
  TVBLimiterEuler(name),
  m_eulerVarSet(CFNULL),
  m_gammaMinusOne()
{
}

//////////////////////////////////////////////////////////////////////////////

TVBLimiterEuler2D::~TVBLimiterEuler2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2D::configure ( Config::ConfigArgs& args )
{
  TVBLimiter::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2D::computePressFromConsVar(const RealVector& consVar, CFreal& press)
{
  const CFreal& rho  = consVar[0];
  const CFreal& rhoU = consVar[1];
  const CFreal& rhoV = consVar[2];
  const CFreal& rhoE = consVar[3];

  press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2D::computeRhoEFromPressAndOtherConsVar(RealVector& consVar, const CFreal& press)
{
  const CFreal& rho  = consVar[0];
  const CFreal& rhoU = consVar[1];
  const CFreal& rhoV = consVar[2];
  CFreal& rhoE = consVar[3];

  rhoE = press/m_gammaMinusOne + 0.5*(rhoU*rhoU+rhoV*rhoV)/rho;
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2D::setup()
{
  CFAUTOTRACE;

  TVBLimiterEuler::setup();

  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in TVBLimiterEuler2DFluxReconstruction!");
  }
  cf_assert(m_nbrEqs == 4);

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiterEuler2D::unsetup()
{
  CFAUTOTRACE;

  TVBLimiterEuler::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
