#include "Chemistry/CH4/CH4.hh"
#include "ChemCH43DPrimArrhenius.hh"
#include "ChemCH43DSourceVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ChemCH43DPrimArrhenius, SourceVarSet, CH4Module, 1> ChemCH43DPrimArrheniusProvider("ChemCH43DPrimArrhenius");

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrimArrhenius::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("m","Exponent coefficient m of the CH4 concentration");
   options.addConfigOption< CFreal >("n","Exponent coefficient n of the O2 concentration");
   options.addConfigOption< CFreal >("A","Arrhenius coefficient");
   options.addConfigOption< CFreal >("EaR","Exponent coefficient of the form Ea/R");
}

//////////////////////////////////////////////////////////////////////////////

ChemCH43DPrimArrhenius::ChemCH43DPrimArrhenius(const std::string& name) :
  ChemCH43DSourceVarSet(name)
{
   addConfigOptionsTo(this);

  _A = 1.3E8;
   setParameter("A",&_A);



  _EaR = 24358.0;
   setParameter("EaR",&_EaR);



  _m = -0.3;
   setParameter("m",&_m);



  _n = 1.3;
   setParameter("n",&_n);



  _multP.resize(getModel()->getNbEquations());
  _multP[XCH4] = 1.0;
  _multP[XO2 ] = 2.0;
  _multP[XCO2] =-1.0;
  _multP[XH2O] =-2.0;

}

//////////////////////////////////////////////////////////////////////////////

ChemCH43DPrimArrhenius::~ChemCH43DPrimArrhenius()
{
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrimArrhenius::configure ( Config::ConfigArgs& args )
{
  ChemCH43DSourceVarSet::configure(args);

  // add here configuration, specific of this class
  const CFreal R = getModel()->getR();;

  _Ea = _EaR * R;
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrimArrhenius::getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
{
  coef = 0.;
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrimArrhenius::getIndepSourceCoefs(const Framework::State& state, RealVector& coef)
{
  // get the temperature
  const CFreal T = getModel()->getTemperature();

  // get the R constant of gases
  const CFreal R = getModel()->getR();;

  // get the pressure
  const CFreal p = getModel()->getPressure();

  const CFreal invRT = 1.0 / (R * T);

  // O2 and CH4 concentration
  const CFreal concCH4 = p * state[XCH4] * invRT;
  const CFreal concO2  = p * state[XO2 ] * invRT;

  const CFreal arrh = - _A * exp( - _Ea * invRT ) * pow(concCH4,_m) * pow(concO2,_n);

  coef = arrh * _multP;
}

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
