#ifndef COOLFluiD_Numerics_FiniteVolume_RMSJouleHeatSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_RMSJouleHeatSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
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
 * This class represents a time-averaged Joule heat source term
 *
 * @author Radek Honzatko
 *
 */
class RMSJouleHeatSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  RMSJouleHeatSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~RMSJouleHeatSourceTerm();

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
   * Set the variable set
   * @pre the input pointer is non const to allow dynamic_cast
   */
  void setVarSet(Common::SafePtr<Framework::ConvectiveVarSet> varSet)
  {
  }

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
  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_rmsJouleHeatSource;

}; // end of class RMSJouleHeatSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RMSJouleHeatSourceTerm_hh
