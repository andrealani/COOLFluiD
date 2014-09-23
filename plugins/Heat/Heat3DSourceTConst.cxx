#include "Heat/Heat.hh"
#include "Heat3DSourceTConst.hh"
#include "Heat3DSourceVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Heat3DSourceTConst, SourceVarSet, HeatModule, 1> heat3DSourceTConstProvider("Heat3DSourceTConst");

//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTConst::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("IndepCoef","Independent Coefficient");
   options.addConfigOption< CFreal >("LinearCoef","Linear Coefficient");
}

//////////////////////////////////////////////////////////////////////////////

Heat3DSourceTConst::Heat3DSourceTConst(const std::string& name) :
  Heat3DSourceVarSet(name)
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

Heat3DSourceTConst::~Heat3DSourceTConst()
{
}

//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTConst::getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
{
  coef = _linearCoef;
}

//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTConst::getIndepSourceCoefs(const Framework::State& state, RealVector& coef)
{
  coef = _indepCoef;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
