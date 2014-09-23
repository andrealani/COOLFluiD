#include "DiffuseReflector.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"


namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<DiffuseReflector,
Reflector,
RadiativeTransferModule,
1>
DiffuseReflectorProvider("DiffuseReflector");

DiffuseReflector::DiffuseReflector(const std::string& name):
    Reflector(name)
{
  addConfigOptionsTo(this);
}

CFreal DiffuseReflector::getReflectionProbability(CFreal lambda, RealVector &s_o, RealVector &s_i, RealVector &normal)
{
 return 0.025330295910584444; //(1./(4.*pi^2)
}

void DiffuseReflector::getRandomDirection(CFreal &lambda, RealVector &s_o, RealVector &s_i, RealVector &normal)
{
  cf_assert(s_o.size() == s_i.size());
  //diffuse: just an emmission
  m_rand.hemiDirections(s_i.size(), normal, s_o );
}




}
}
