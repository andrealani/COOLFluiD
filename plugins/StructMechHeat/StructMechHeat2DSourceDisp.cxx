#include "StructMechHeat/StructMechHeat.hh"
#include "StructMechHeat2DSourceDisp.hh"
#include "StructMechHeat2DSourceVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMechHeat2DSourceDisp, SourceVarSet, StructMechHeatModule, 1> StructMechHeat2DSourceDispProvider("StructMechHeat2DSourceDisp");

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DSourceDisp::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("LinearDef","Definition of the Linear Coefficient Functions.");
   options.addConfigOption< std::vector<std::string> >("IndepDef","Definition of the Independent Coefficient Functions.");
}

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2DSourceDisp::StructMechHeat2DSourceDisp(const std::string& name) :
  StructMechHeat2DSourceVarSet(name)
{
   addConfigOptionsTo(this);
  vector<std::string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "T";
  setVarNames(names);

  _functionsIndep = std::vector<std::string>();
   setParameter("IndepDef",&_functionsIndep);

  _functionsLinear = std::vector<std::string>();
   setParameter("LinearDef",&_functionsLinear);

  _varsDescript.resize(6);
  _varsDescript[0] = "x";
  _varsDescript[1] = "y";
  _varsDescript[2] = "rho";
  _varsDescript[3] = _varNames[0];
  _varsDescript[4] = _varNames[1];
  _varsDescript[5] = _varNames[2];

  _vars.resize(_varsDescript.size());

}

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2DSourceDisp::~StructMechHeat2DSourceDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DSourceDisp::configure ( Config::ConfigArgs& args )
{
  StructMechHeat2DSourceVarSet::configure(args);

  // add here configuration, specific of this class

  if (hasIndepCoef()){
///@todo this cannot be!!!
  // configure the Independent Function
//  cf_assert(_functionsIndep.size() == PhysicalModelStack::getActive()->getNbEq());
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
///@todo this cannot be!!!
  // configure the Linear Function
//  cf_assert(_functionsLinear.size() == PhysicalModelStack::getActive()->getNbEq());
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

void StructMechHeat2DSourceDisp::setup()
{
   CFAUTOTRACE;

   StructMechHeat2DSourceVarSet::setup();

  _vectorialCoef.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DSourceDisp::getLinearSourceCoefs(const State& state, RealMatrix& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[2]  = getModel()->getDensity();
  _vars[3]  = state[0];
  _vars[4]  = state[1];
  _vars[5]  = state[2];

  ///@todo CHANGE THIS!!
///@todo This is not correct
  if (!_functionsLinear.empty()){
    _vFunctionLinear.evaluate(_vars,_vectorialCoef);
  }

cf_assert(false);
}


//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DSourceDisp::getIndepSourceCoefs(const State& state, RealVector& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[2]  = getModel()->getDensity();
  _vars[3]  = state[0];
  _vars[4]  = state[1];
  _vars[5]  = state[2];

 ///@todo CHANGE THIS!!
  if (!_functionsIndep.empty()){
   _vFunctionIndep.evaluate(_vars,coef);
   }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
