#ifndef COOLFluiD_Numerics_FiniteVolume_SuperOutletMHD3DProjectionHeating_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperOutletMHD3DProjectionHeating_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic outlet command in 3D
   * for projection scheme on volumetric heating modelling of the solar wind
   * 
   * @author Andrea Lani
   * @author Mehmet Sarp Yalim
   *
   */
class SuperOutletMHD3DProjectionHeating : public FVMCC_BC {
public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SuperOutletMHD3DProjectionHeating(const std::string& name);
  
  /**
   * Default destructor
   */
  ~SuperOutletMHD3DProjectionHeating();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object by complementing the 
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: //data

  /// phi value that is to be fixed
  CFreal _refPhi;

  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// socket for the PFSS magnetic field components in Cartesian coordinates to be computed once in the setup phase
  Framework::DataSocketSink<std::vector<CFreal> > socket_BPFSS;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

}; // end of class SuperOutletMHD3DProjectionHeating

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperOutletMHD3DProjectionHeating_hh
