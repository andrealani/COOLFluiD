#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "TwoLayerWeakFarFieldEuler2DConsAlphaMTP.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerWeakFarFieldEuler2DConsAlphaMTP, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> TwoLayerWeakFarFieldEuler2DConsAlphaMTPProvider("TwoLayerWeakFarFieldEuler2DConsAlphaMTP");

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarFieldEuler2DConsAlphaMTP::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("mach","Mach number");
   options.addConfigOption< CFreal >("alpha","flow angle");
   options.addConfigOption< CFreal >("T","temperature");
   options.addConfigOption< CFreal >("p","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakFarFieldEuler2DConsAlphaMTP::TwoLayerWeakFarFieldEuler2DConsAlphaMTP(const std::string& name) :
  TwoLayerWeakBC2D(name),
  _varSet()
 {
   addConfigOptionsTo(this);
   _alpha = 0.0;
   setParameter("alpha",&_alpha);

   _mach = 0.0;
   setParameter("mach",&_mach);

   _temperature = 0.0;
   setParameter("T",&_temperature);

  _pressure = 0.0;
   setParameter("p",&_pressure);
 }

//////////////////////////////////////////////////////////////////////////////

TwoLayerWeakFarFieldEuler2DConsAlphaMTP::~TwoLayerWeakFarFieldEuler2DConsAlphaMTP()
{
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarFieldEuler2DConsAlphaMTP::setup()
{
  TwoLayerWeakBC2D::setup();
  _varSet->setup();
  // convert angle in radiants (/// @todo do the same in FVM BCs)
  _alpha = _alpha*MathTools::MathConsts::CFrealPi()/180.;

  _temperature /= _varSet->getModel()->getTempRef();
  _pressure /= _varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarFieldEuler2DConsAlphaMTP::configure ( Config::ConfigArgs& args )
{
  TwoLayerWeakBC2D::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerWeakFarFieldEuler2DConsAlphaMTP::setGhostState(const State& state,
                                                    State& gstate)
{
  const CFreal R = _varSet->getModel()->getR();
  const CFreal rho = _pressure/(R*_temperature);
  const CFreal tgAlpha = tan(_alpha);
  const CFreal u = _mach*sqrt(_varSet->getModel()->getGamma()*R*
                              _temperature/(1. + tgAlpha*tgAlpha));

  gstate[0] = rho;
  gstate[1] = rho*u;
  gstate[2] = gstate[1]*tgAlpha;

  const CFreal rhoV2 = gstate[1]*gstate[1] + gstate[2]*gstate[2];
  gstate[3] = _pressure/(_varSet->getModel()->getGamma() - 1.) +
    0.5*rhoV2/rho;
 }

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
