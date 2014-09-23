#include "Heat/Heat.hh"
#include "Heat3DSourceTDep.hh"
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

Environment::ObjectProvider<Heat3DSourceTDep, SourceVarSet, HeatModule, 1> heat3DSourceTDepProvider("Heat3DSourceTDep");

//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTDep::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("LinearDef","Definition of the Linear Coefficient Functions.");
   options.addConfigOption< std::vector<std::string> >("IndepDef","Definition of the Independent Coefficient Functions.");
}

//////////////////////////////////////////////////////////////////////////////

Heat3DSourceTDep::Heat3DSourceTDep(const std::string& name) :
  Heat3DSourceVarSet(name)
{

  addConfigOptionsTo(this);
  vector<std::string> names(1);
  names[0] = "T";
  setVarNames(names);

  _functionsIndep = std::vector<std::string>();
  setParameter("IndepDef",&_functionsIndep);

  _functionsLinear = std::vector<std::string>();
  setParameter("LinearDef",&_functionsLinear);

  _varsDescript.resize(4);
  _varsDescript[0] = "x";
  _varsDescript[1] = "y";
  _varsDescript[2] = "z";
  _varsDescript[3] = _varNames[0];

  _vars.resize(_varsDescript.size());

}

//////////////////////////////////////////////////////////////////////////////

Heat3DSourceTDep::~Heat3DSourceTDep()
{
}

//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTDep::setup()
{
   Heat3DSourceVarSet::setup();

  _vectorialCoef.resize(PhysicalModelStack::getActive()->getNbEq());
}


//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTDep::configure ( Config::ConfigArgs& args )
{
  Heat3DSourceVarSet::configure(args);

  // add here configuration, specific of this class

  // configure the Independent Function
  if (hasIndepCoef()) {
///@todo this cannot be!!!
//    cf_assert(_functionsIndep.size() == PhysicalModelStack::getActive()->getNbEq());
    _vFunctionIndep.setFunctions(_functionsIndep);
    _vFunctionIndep.setVariables(_varsDescript);
    try {
      _vFunctionIndep.parse();
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
  }

  // configure the Linear Function
  if (hasLinearCoef()) {
///@todo this cannot be!!!
//    cf_assert(_functionsLinear.size() == PhysicalModelStack::getActive()->getNbEq());
    _vFunctionLinear.setFunctions(_functionsLinear);
    _vFunctionLinear.setVariables(_varsDescript);
    try {
      _vFunctionLinear.parse();
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTDep::getLinearSourceCoefs(const State& state, RealMatrix& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[ZZ] = coord[ZZ];
  _vars[3]  = state[0];

  _vFunctionLinear.evaluate(_vars,_vectorialCoef);
  coef(0,0)=_vectorialCoef[0];
}


//////////////////////////////////////////////////////////////////////////////

void Heat3DSourceTDep::getIndepSourceCoefs(const State& state, RealVector& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[ZZ] = coord[ZZ];
  _vars[3]  = state[0];
  _vFunctionIndep.evaluate(_vars,coef);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
