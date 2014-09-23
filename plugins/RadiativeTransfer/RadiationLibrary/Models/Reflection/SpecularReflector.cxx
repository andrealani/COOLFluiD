#include "SpecularReflector.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathFunctions.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<SpecularReflector,
                            Reflector,
                            RadiativeTransferModule,
                            1>
SpecularReflectorProvider("SpecularReflector");

SpecularReflector::SpecularReflector(const std::string& name):
     Reflector(name)
{
  addConfigOptionsTo(this);
}

CFreal SpecularReflector::getReflectionProbability(CFreal lambda, RealVector &s_o, RealVector &s_i, RealVector &normal)
{
  cf_assert(s_o.size() == s_i.size());
  return (s_i==s_o) ? 1. : 0.;
}

void SpecularReflector::getRandomDirection(CFreal &lambda, RealVector &s_o, RealVector &s_i, RealVector &normal)
{
  cf_assert(s_o.size() == s_i.size());
  cf_assert(s_o.size() == normal.size());

  //specular: snell law
  s_o=s_i-2.0*MathTools::MathFunctions::innerProd(normal,s_i)*normal;
}


}
}
