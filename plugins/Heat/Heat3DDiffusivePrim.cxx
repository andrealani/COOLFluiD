#include "Environment/ObjectProvider.hh"

#include "Heat/Heat.hh"
#include "Heat/Heat3DDiffusivePrim.hh"
#include "Heat/Heat3DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Heat3DDiffusivePrim,
                            DiffusiveVarSet,
                            HeatModule, 2>
aHeat3DDiffusivePrimProvider("Heat3DDiffusivePrim");

//////////////////////////////////////////////////////////////////////////////

Heat3DDiffusivePrim::Heat3DDiffusivePrim(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Heat3DDiffusiveVarSet(name, model)
{
  vector<std::string> names(1);
  names[0] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Heat3DDiffusivePrim::~Heat3DDiffusivePrim()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
