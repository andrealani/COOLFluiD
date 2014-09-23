#include "PostProcessNull.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<PostProcessNull, PostProcess, RadiativeTransferModule, 1 >
PostProcessNullProvider("PostProcessNull");


}
}
