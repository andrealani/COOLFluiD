#include "StructMech/StructMech.hh"
#include "StructMech3DInertiaDisp.hh"
#include "StructMech3DInertiaVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech3DInertiaDisp, InertiaVarSet, StructMechModule, 1>
structMech3DInertiaDispProvider("StructMech3DInertiaDisp");

//////////////////////////////////////////////////////////////////////////////

StructMech3DInertiaDisp::StructMech3DInertiaDisp(const std::string& name) :
  StructMech3DInertiaVarSet(name)
{
  vector<std::string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "w";
  setVarNames(names);

}

//////////////////////////////////////////////////////////////////////////////

StructMech3DInertiaDisp::~StructMech3DInertiaDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech3DInertiaDisp::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  StructMech3DInertiaVarSet::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void StructMech3DInertiaDisp::setup()
{
  StructMech3DInertiaVarSet::setup();

  _rho = getModel()->getDensity();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
