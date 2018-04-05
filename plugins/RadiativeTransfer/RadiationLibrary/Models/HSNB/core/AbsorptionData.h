#ifndef COOLFluiD_RadiativeTransfer_ABSORPTIONDATA_H
#define COOLFluiD_RadiativeTransfer_ABSORPTIONDATA_H

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"


namespace COOLFluiD {

    namespace RadiativeTransfer {


    enum trajectoryState {NOTASSIGNED=-1};
//////////////////////////////////////////////////////////////////////////////
/// \brief The PhotonData struct
///


struct AbsorptionData{

    CFreal traceDistance;
    
    
};

//////////////////////////////////////////////////////////////////////////////

    }
}

#endif // ABSORPTIONDATA_H
