#ifndef COOLFluiD_RadiativeTransfer_PhotonData_hh
#define COOLFluiD_RadiativeTransfer_PhotonData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"


namespace COOLFluiD {

    namespace RadiativeTransfer {
      
    enum trajectoryState {NOTASSIGNED=-1};
      
//////////////////////////////////////////////////////////////////////////////

/// \brief The PhotonData struct
///
      
struct PhotonData {
  
  CFint globalTrajectoryId=NOTASSIGNED;
  
  /**
   * @brief used to determine the ordering of trajectories
   *
   * Will be 1 for the partition of the process that spawned the photon and increased
     * every time the photon passes overt to a new process.
     */
  CFint indexWithinTrajectory=0;
  
  /**
   * @brief rank of the proc. that spawned this photon.
   */
  CFint fatherProcessId;
  CFreal KS;
  CFreal energyFraction;
  CFreal wavelength;
};
      
//////////////////////////////////////////////////////////////////////////////
      
    }
}

#endif // PHOTONDATA_H
