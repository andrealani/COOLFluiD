#ifndef COOLFluiD_Physics_Maxwell_Maxwell3DProjectionConsT_hh
#define COOLFluiD_Physics_Maxwell_Maxwell3DProjectionConsT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Maxwell3DProjectionVarSetT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a GPU-enabled Maxwell physical model 2D for projection scheme
 * for conservative variables
 *
 * @author Alejandro Alvarez Laguna
 * @author Andrea Lani
 * @author Isaac Alonso
 */
class Maxwell3DProjectionConsT : public Maxwell3DProjectionVarSetT {
public: 

  /**
   * Constructor
   */
  HOST_DEVICE Maxwell3DProjectionConsT(MaxwellProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    Maxwell3DProjectionVarSetT(dco) {}
 
   /**
   * Constructor
   */
   HOST_DEVICE Maxwell3DProjectionConstT() : Maxwell3DProjectionVarSetT() {}  



   /**
   * Default destructor
   */
  HOST_DEVICE ~Maxwell3DProjectionConsT();

  

  HOST_DEVICE void computePhysicalData(CFreal* state, CFreal* data)
  {  
    //CFLog(NOTICE, "computePhysicalData\n");  
    data[ConvMaxwellTerm::BX] = state[0];
    data[ConvMaxwellTerm::BY] = state[1];
    data[ConvMaxwellTerm::BZ] = state[2];
    data[ConvMaxwellTerm::EX] = state[3];
    data[ConvMaxwellTerm::EY] = state[4];
    data[ConvMaxwellTerm::EZ] = state[5];
    data[MaxwellProjectionTerm::PSI] = state[6];  
    data[MaxwellProjectionTerm::PHI] = state[7];
  }

//////////////////////////////////////////////////////////////////////////////

  HOST_DEVICE void computeStateFromPhysicalData(CFreal* data, CFreal* state)
  {
    state[0] = data[ConvMaxwellTerm::BX];
    state[1] = data[ConvMaxwellTerm::BY];
    state[2] = data[ConvMaxwellTerm::BZ];
    state[3] = data[ConvMaxwellTerm::EX];
    state[4] = data[ConvMaxwellTerm::EY];
    state[5] = data[ConvMaxwellTerm::EZ];
    state[6] = data[MaxwellProjectionTerm::PSI];
    state[7] = data[MaxwellProjectionTerm::PHI];
  }
 

 
}; // end of class Maxwell3DProjectionConsT

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_Maxwell3DProjectionConsT_hh
