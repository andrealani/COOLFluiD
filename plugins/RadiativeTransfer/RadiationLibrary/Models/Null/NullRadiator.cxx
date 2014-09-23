#include "NullRadiator.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<NullRadiator,
                            Radiator,
                            RadiativeTransferModule,
                            1>
NullRadiatorProvider("NullRadiator");

}
}
