#ifndef COOLFluiD_RadiativeTransfer_HSNBPhotonData_hh
#define COOLFluiD_RadiativeTransfer_HSNBPhotonData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include <vector>
#include "RadiativeTransfer/Solvers/MonteCarlo/PhotonData.hh"

namespace COOLFluiD {

    namespace RadiativeTransfer {


//////////////////////////////////////////////////////////////////////////////
/// \brief The PhotonData struct
///
const int HSNBBUFSIZE  = 1000;

enum HSNBMechanismType {THICKDIATOMIC=0, THINDIATOMIC=1, CONTINUUM=2,
                        ATOMICLINE=3, BASETYPE=4, WALLEMISSION=5};

struct HSNBPhotonData: PhotonData {

    CFuint targetProcessId;
    CFuint prevProcessId;
    CFint mechType;
    //NEEDED FOR ALL MECHANISMS
    CFuint mechanismIndex;
    CFuint nbCrossedCells;
    CFreal energyResiduum;
    CFreal transm;
    CFreal wallTransmissivity;


};




    }
}

#endif // PHOTONDATA_H
