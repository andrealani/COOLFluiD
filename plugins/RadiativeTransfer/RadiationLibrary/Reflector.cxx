#include "Reflector.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

/// Constructor without arguments
Reflector::Reflector(const std::string& name):
           Common::OwnedObject(),
           ConfigObject(name)
{
  addConfigOptionsTo(this);
}


}
}
