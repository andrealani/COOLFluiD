#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeSourceTermFVMCC_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeSourceTermFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeSourceTerm.hh"
#include "Framework/Storage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
    class ConvectiveVarSet;
    class DiffusiveVarSet;
  }

  namespace Numerics {

    namespace FiniteVolume {
      class CellCenterFVMData;
      
//////////////////////////////////////////////////////////////////////////////

 /**
  * This class offers and abstract interface to computes the source
  * term for cell center FVM schemes.
  *
  * @author Andrea Lani
  *
  */
class ComputeSourceTermFVMCC : public Framework::ComputeSourceTerm<CellCenterFVMData> {

public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  ComputeSourceTermFVMCC(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ComputeSourceTermFVMCC();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian) = 0;

  /// Returns the DataSocket's that this strategy provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /// get the velocity IDs in the State (update variables)
  const std::vector<CFuint>& getStateVelocityIDs() const {return m_velIDs;}
  
protected:

  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// socket for the normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// socket for the states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// socket for the isOutward flag
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// socket for the gradient ux
  Framework::DataSocketSink<CFreal> socket_uX;
  
  /// socket for the gradient uy
  Framework::DataSocketSink<CFreal> socket_uY;
  
  /// socket for the gradient uz
  Framework::DataSocketSink<CFreal> socket_uZ;
  
  /// data handle cashed for efficiency  reasons
  Framework::DataHandle<CFreal> m_ux;
  
  /// data handle cashed for efficiency  reasons
  Framework::DataHandle<CFreal> m_uy;
  
  /// data handle cashed for efficiency  reasons
  Framework::DataHandle<CFreal> m_uz;
  
  /// gradients exist
  bool m_gradientsExist;
  
  /// IDs corresponding to the velocity components
  std::vector<CFuint> m_velIDs;
  
  /// physical data array
  RealVector m_pdataArray;
  
  /// use the gradient computed with the least square reconstruction
  bool m_useGradientLS;
  
}; // end of class ComputeSourceTermFVMCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeSourceTermFVMCC_hh
