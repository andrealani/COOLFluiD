#include "FEMIntegrableEntity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

FEMIntegrableEntity::FEMIntegrableEntity(const std::string& name) : Framework::MethodStrategy<FiniteElementMethodData>(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FEMIntegrableEntity::~FEMIntegrableEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void FEMIntegrableEntity::setup()
{
  FiniteElementMethodStrategy::setup();

  // add specific configuration here
  _localElemData = &(getMethodData().getLocalElementData());
}

//////////////////////////////////////////////////////////////////////////////

void FEMIntegrableEntity::integrateContourInGeo(GeometricEntity* const cell,
                                             RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"FEMIntegrableEntity::integrateContourInGeo()");
}

//////////////////////////////////////////////////////////////////////////////

void FEMIntegrableEntity::integrateVolumeInGeo(GeometricEntity* const cell,
                                            RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"FEMIntegrableEntity::integrateVolumeInGeo()");
}

//////////////////////////////////////////////////////////////////////////////

void FEMIntegrableEntity::integrateContourInGeo(GeometricEntity* const cell,
                                             RealMatrix& result)
{
  throw Common::NotImplementedException (FromHere(),"FEMIntegrableEntity::integrateContourInGeo()");
}

//////////////////////////////////////////////////////////////////////////////

void FEMIntegrableEntity::integrateVolumeInGeo(GeometricEntity* const cell,
                                            RealMatrix& result)
{
  throw Common::NotImplementedException (FromHere(),"FEMIntegrableEntity::integrateVolumeInGeo()");
}

//////////////////////////////////////////////////////////////////////////////
    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
