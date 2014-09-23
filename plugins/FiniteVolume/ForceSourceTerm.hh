#ifndef COOLFluiD_Numerics_FiniteVolume_ForceSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_ForceSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a momentum source term ... still under devt.
 *
 * @author Frederic Vandermot
 *
 */
class ForceSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  ForceSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ForceSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
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
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data
  ///// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_ForceSourceXcomponent;
  Framework::DataSocketSink<CFreal> socket_ForceSourceYcomponent;
  Framework::DataSocketSink<CFreal> socket_ForceSourceZcomponent;

}; // end of class ForceSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ForceSourceTerm_hh
