#ifndef COOLFluiD_Numerics_FiniteVolume_DiffFluxCarbuncleFixEuler_hh
#define COOLFluiD_Numerics_FiniteVolume_DiffFluxCarbuncleFixEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeDiffusiveFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model
 *
 * @author Andrea Lani
 *
 */
class DiffFluxCarbuncleFixEuler : public ComputeDiffusiveFlux {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  DiffFluxCarbuncleFixEuler(const std::string& name);

  /**
   * Default destructor
   */
  ~DiffFluxCarbuncleFixEuler();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up 
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result);

protected:
  
  /**
   * Set the quantities needed to compute gradients (velocity magnitude, temperature, Mach number)
   * starting from the states
   */
  void setGradientVars(const std::vector<RealVector*>& states,
		       RealMatrix& values,
		       CFuint stateSize);
  
protected: // data
  
  /// socket for the cell volumes storage
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// socket indicating where the carbuncle fix is active
  Framework::DataSocketSource<CFreal> socket_fix_active;
  
  /// number states in the control volume
  CFuint _nbCVStates;
  
  // array of states involved in the flux computation
  std::vector<RealVector*> _states;
  
  // array of values (p, u, v, T, ...)
  RealMatrix _values;
  
  /// arrray of gradients
  std::vector<RealVector*> _gradients;
  
  // array of average state in update variables
  RealVector _avState;
  
  // array of average physical data in update variables
  RealVector _avPhysData;
  
  /// articficial viscosity
  CFreal _muS;
  
  /// dissipation coefficient
  CFreal _eps; 
  
}; // end of class DiffFluxCarbuncleFixEuler
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DiffFluxCarbuncleFixEuler_hh
