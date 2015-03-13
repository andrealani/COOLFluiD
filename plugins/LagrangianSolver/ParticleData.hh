#include "Common/COOLFluiD.hh"
#include "RadiativeTransfer/RadiativeTransferModule.hh"

/**
 * This structs model a Particle
 *
 * @author Pedro Santos
 * @author Alessandro Sanna
 * @author Andrea Lani
 */
namespace COOLFluiD {

namespace LagrangianSolver {

struct CommonData{
    CFreal direction[3];
    CFreal currentPoint[3];
    CFuint cellID;
};


template<class UserData1>
struct Particle
{
  CommonData commonData;
  UserData1 userData;
};
}

}
