#ifndef COOLFluiD_Physics_Maxwell_Maxwell2DProjectionConsT_hh
#define COOLFluiD_Physics_Maxwell_Maxwell2DProjectionConsT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Maxwell2DProjectionVarSetT.hh"

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
class Maxwell2DProjectionConsT : public Maxwell2DProjectionVarSetT {
public: // function

   /**
   * Constructor
   */
  HOST_DEVICE Maxwell2DProjectionConsT(MaxwellProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    Maxwell2DProjectionVarSetT(dco) {}
 
  /**
   * Constructor
   */
  HOST_DEVICE Maxwell2DProjectionConsT() : Maxwell2DProjectionVarSetT() {}



  /**
   * Default destructor
   */
  HOST_DEVICE ~Maxwell2DProjectionConsT() {};

  

  HOST_DEVICE void computePhysicalData(CFreal* state, CFreal* data)
  {  
    
    //CFLog(VERBOSE, "computePhysicalData\n");  
    data[0] = state[0];
    data[1] = state[1];
    data[2] = state[2];
    data[3] = state[3];
    data[4] = state[4];
    data[5] = state[5];
    data[6] = state[6];  
    data[7] = state[7];
  }

//////////////////////////////////////////////////////////////////////////////

  HOST_DEVICE void computeStateFromPhysicalData(CFreal* data, CFreal* state) 
  {
    state[0] = data[0];
    state[1] = data[1];
    state[2] = data[2];
    state[3] = data[3];
    state[4] = data[4];
    state[5] = data[5];
    state[6] = data[6];
    state[7] = data[7];
  }
 

 
}; // end of class Maxwell2DProjectionConsT

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_Maxwell2DProjectionConsT_hh
