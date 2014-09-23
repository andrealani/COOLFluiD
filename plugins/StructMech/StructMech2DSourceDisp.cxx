#include "StructMech/StructMech.hh"
#include "StructMech2DSourceDisp.hh"
#include "StructMech2DSourceVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech2DSourceDisp, SourceVarSet, StructMechModule, 1> structMech2DSourceDispProvider("StructMech2DSourceDisp");

//////////////////////////////////////////////////////////////////////////////

void StructMech2DSourceDisp::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("LinearDef","Definition of the Linear Coefficient Functions.");
   options.addConfigOption< std::vector<std::string> >("IndepDef","Definition of the Independent Coefficient Functions.");
}

//////////////////////////////////////////////////////////////////////////////

StructMech2DSourceDisp::StructMech2DSourceDisp(const std::string& name) :
  StructMech2DSourceVarSet(name)
{
   addConfigOptionsTo(this);
  vector<std::string> names(2);
  names[0] = "u";
  names[1] = "v";
  setVarNames(names);

  _functionsIndep = std::vector<std::string>();
   setParameter("IndepDef",&_functionsIndep);

  _functionsLinear = std::vector<std::string>();
   setParameter("LinearDef",&_functionsLinear);

  _varsDescript.resize(5);
  _varsDescript[0] = "x";
  _varsDescript[1] = "y";
  _varsDescript[2] = "rho";
  _varsDescript[3] = _varNames[0];
  _varsDescript[4] = _varNames[1];

  _vars.resize(_varsDescript.size());

}

//////////////////////////////////////////////////////////////////////////////

StructMech2DSourceDisp::~StructMech2DSourceDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech2DSourceDisp::configure ( Config::ConfigArgs& args )
{
  StructMech2DSourceVarSet::configure(args);

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

void StructMech2DSourceDisp::setup()
{
   CFAUTOTRACE;

   StructMech2DSourceVarSet::setup();

  _vectorialCoef.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void StructMech2DSourceDisp::getLinearSourceCoefs(const State& state, RealMatrix& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[2]  = getModel()->getDensity();
  _vars[3]  = state[0];
  _vars[4]  = state[1];

  ///@todo CHANGE THIS!!
  if (!_functionsLinear.empty()){
    _vFunctionLinear.evaluate(_vars,_vectorialCoef);
    coef(0,0)=_vectorialCoef[0];
    }
}


//////////////////////////////////////////////////////////////////////////////

void StructMech2DSourceDisp::getIndepSourceCoefs(const State& state, RealVector& coef)
{
  RealVector& coord = state.getCoordinates();

  _vars[XX] = coord[XX];
  _vars[YY] = coord[YY];
  _vars[2]  = getModel()->getDensity();
  _vars[3]  = state[0];
  _vars[4]  = state[1];

 ///@todo CHANGE THIS!!
 if (!_functionsIndep.empty()){
  _vFunctionIndep.evaluate(_vars,coef);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
