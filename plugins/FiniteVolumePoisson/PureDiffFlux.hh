#ifndef COOLFluiD_Numerics_FiniteVolume_PureDiffFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_PureDiffFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeDiffusiveFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace Poisson {
      class PoissonDiffVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Poisson model
 *
 * @author Alejandro Alvarez Laguna
 *
 */
class PureDiffFlux : public ComputeDiffusiveFlux {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  PureDiffFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~PureDiffFlux();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up 
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result);

protected: // data
  
  /// socket for the cell volumes storage
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// number states in the control volume
  CFuint _nbCVStates;
  
  // array of states involved in the flux computation
  std::vector<RealVector*> _states;

  // array of values (phi)
  RealMatrix _values;
  
  /// arrray of gradients
  std::vector<RealVector*> _gradients;
  
  // array of average state in update variables
  RealVector _avState;

  /// VarSet
  Common::SafePtr<Physics::Poisson::PoissonDiffVarSet> _varSet;
    
}; // end of class PureDiffFlux
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PureDiffFlux_hh
