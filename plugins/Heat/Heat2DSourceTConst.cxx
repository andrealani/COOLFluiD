#include "Heat/Heat.hh"
#include "Heat2DSourceTConst.hh"
#include "Heat2DSourceVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Heat2DSourceTConst, SourceVarSet, HeatModule, 1> heat2DSourceTConstProvider("Heat2DSourceTConst");

//////////////////////////////////////////////////////////////////////////////

void Heat2DSourceTConst::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("IndepCoef","Independent Coefficient");
   options.addConfigOption< CFreal >("LinearCoef","Linear Coefficient");
}

//////////////////////////////////////////////////////////////////////////////

Heat2DSourceTConst::Heat2DSourceTConst(const std::string& name) :
  Heat2DSourceVarSet(name)
{
   addConfigOptionsTo(this);
  vector<std::string> names(1);
  names[0] = "T";
  setVarNames(names);

  _linearCoef = 0.;
   setParameter("LinearCoef",&_linearCoef);



  _indepCoef = 0.;
   setParameter("IndepCoef",&_indepCoef);

}

//////////////////////////////////////////////////////////////////////////////

Heat2DSourceTConst::~Heat2DSourceTConst()
{
}

//////////////////////////////////////////////////////////////////////////////

void Heat2DSourceTConst::getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
{
  coef = _linearCoef;
}

//////////////////////////////////////////////////////////////////////////////

void Heat2DSourceTConst::getIndepSourceCoefs(const Framework::State& state, RealVector& coef)
{
  coef = _indepCoef;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
