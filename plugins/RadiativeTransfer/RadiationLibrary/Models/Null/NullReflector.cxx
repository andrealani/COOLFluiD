#include "NullReflector.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"


namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<NullReflector,
                            Reflector,
                            RadiativeTransferModule,
                           1>
NullReflectorProvider("NullReflector");

NullReflector::NullReflector(const std::string& name):
    Reflector(name)
{
  addConfigOptionsTo(this);
}




}
}
