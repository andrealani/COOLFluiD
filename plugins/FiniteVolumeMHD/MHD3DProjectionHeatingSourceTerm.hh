#ifndef COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionHeatingSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionHeatingSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an additional volumetric heating source term for a 
 * more improved  solar wind model than the polytropic model implemented 
 * from Francesco Zuccarello's PhD thesis for 3D conservative variables used 
 * with the hyperbolic divergence cleaning method
 *
 * @author Mehmet Sarp Yalim
 *
 */
class MHD3DProjectionHeatingSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  MHD3DProjectionHeatingSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD3DProjectionHeatingSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
    
    _globalSockets.createSocketSink<Framework::State*>("states");
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
 
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
  
private: // data

  /// corresponding variable set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// MHD physical data
  RealVector _physicalData;

  /// MHD physical data
  RealVector _dataLeftState;

   /// MHD physical data
  RealVector _dataRightState;

}; // end of class MHD3DProjectionHeatingSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionHeatingSourceTerm_hh
