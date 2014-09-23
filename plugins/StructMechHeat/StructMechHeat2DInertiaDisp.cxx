#include "StructMechHeat/StructMechHeat.hh"
#include "StructMechHeat2DInertiaDisp.hh"
#include "StructMechHeat2DInertiaVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMechHeat2DInertiaDisp, InertiaVarSet, StructMechHeatModule, 1>
StructMechHeat2DInertiaDispProvider("StructMechHeat2DInertiaDisp");

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2DInertiaDisp::StructMechHeat2DInertiaDisp(const std::string& name) :
  StructMechHeat2DInertiaVarSet(name)
{
  vector<std::string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2DInertiaDisp::~StructMechHeat2DInertiaDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DInertiaDisp::setup()
{
  CFAUTOTRACE;

  StructMechHeat2DInertiaVarSet::setup();

  _rho = getModel()->getDensity();
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DInertiaDisp::configure ( Config::ConfigArgs& args )
{
  StructMechHeat2DInertiaVarSet::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
