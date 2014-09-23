#include "StructMech/StructMech.hh"
#include "StructMech3DSourceDisp.hh"
#include "StructMech3DSourceVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech3DSourceDisp, SourceVarSet, StructMechModule, 1> structMech3DSourceDispProvider("StructMech3DSourceDisp");

//////////////////////////////////////////////////////////////////////////////

void StructMech3DSourceDisp::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("LinearDef","Definition of the Linear Coefficient Functions.");
   options.addConfigOption< std::vector<std::string> >("IndepDef","Definition of the Independent Coefficient Functions.");
}

//////////////////////////////////////////////////////////////////////////////

StructMech3DSourceDisp::StructMech3DSourceDisp(const std::string& name) :
  StructMech3DSourceVarSet(name)
{
   addConfigOptionsTo(this);
  vector<std::string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "w";
  setVarNames(names);

  _functionsIndep = std::vector<std::string>();
   setParameter("IndepDef",&_functionsIndep);

  _functionsLinear = std::vector<std::string>();
   setParameter("LinearDef",&_functionsLinear);

  _varsDescript.resize(7);
  _varsDescript[0] = "x";
  _varsDescript[1] = "y";
  _varsDescript[2] = "z";
  _varsDescript[3] = "rho";
  _varsDescript[4] = _varNames[0];
  _varsDescript[5] = _varNames[1];
  _varsDescript[6] = _varNames[2];

  _vars.resize(_varsDescript.size());
  _vectorialCoef.resize(PhysicalModelStack::getActive()->getNbEq());

}

//////////////////////////////////////////////////////////////////////////////

StructMech3DSourceDisp::~StructMech3DSourceDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech3DSourceDisp::configure ( Config::ConfigArgs& args )
{
  StructMech3DSourceVarSet::configure(args);

  // add here configuration, specific of this class

  if (hasIndepCoef()){
  // configure the Independent Function
///@todo this cannot be!!!
  cf_assert(_functionsIndep.size() == PhysicalModelStack::getActive()->getNbEq());
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

  if (hasLinearCoef()){
  // configure the Linear Function
///@todo this cannot be!!!
  cf_assert(_functionsLinear.size() == PhysicalModelStack::getActive()->getNbEq());
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

void StructMech3DSourceDisp::getLinearSourceCoefs(const State& state, RealMatrix& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[ZZ] = coord[ZZ];
  _vars[3]  = getModel()->getDensity();
  _vars[4]  = state[0];
  _vars[5]  = state[1];
  _vars[6]  = state[2];

  /// @todo CHANGE THIS!!
  if (!_functionsLinear.empty()){
    _vFunctionLinear.evaluate(_vars,_vectorialCoef);
    coef(0,0)=_vectorialCoef[0];
  }
}


//////////////////////////////////////////////////////////////////////////////

void StructMech3DSourceDisp::getIndepSourceCoefs(const State& state, RealVector& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[ZZ] = coord[ZZ];
  _vars[3]  = getModel()->getDensity();
  _vars[4]  = state[0];
  _vars[5]  = state[1];
  _vars[6]  = state[2];

  if (!_functionsIndep.empty()){
    _vFunctionIndep.evaluate(_vars,coef);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
