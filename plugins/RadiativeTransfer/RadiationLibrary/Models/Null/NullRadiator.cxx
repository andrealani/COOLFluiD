#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"
#include "RadiativeTransfer/RadiationLibrary/Models/Null/NullRadiator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullRadiator,
                            Radiator,
                            RadiativeTransferModule,
                            1>
NullRadiatorProvider("NullRadiator");

//////////////////////////////////////////////////////////////////////////////

}
}
