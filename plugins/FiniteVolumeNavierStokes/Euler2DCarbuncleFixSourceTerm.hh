#ifndef COOLFluiD_Numerics_FiniteVolume_Euler2DCarbuncleFixSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_Euler2DCarbuncleFixSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DCarbuncleFixSourceTerm : public ComputeSourceTermFVMCC {

public:
  
  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  Euler2DCarbuncleFixSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  ~Euler2DCarbuncleFixSourceTerm();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
    _sockets.createSocketSink<RealVector>("nstates");
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

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
protected:
  
  /**
   * Set the quantities needed to compute gradients (velocity magnitude, temperature, Mach number)
   * starting from the states
   */
  void setGradientVars(const std::vector<RealVector*>& states,
		       RealMatrix& values,
		       CFuint stateSize);

private: // data
  
  /// socket indicating where the carbuncle fix is active
  Framework::DataSocketSource<CFreal> socket_fix_active;
  
  /// socket storing the mu_s
  Framework::DataSocketSource<CFreal> socket_ArtViscCoeff;
  
  /// socket for the normals
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// corresponding variable set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;

  /// vector to store temporary result
  RealVector _temp;

  /// Euler physical data
  RealVector _physicalData;

  /// array for temporary dX vector
  RealVector _dX;

  /// array for temporary speeds
  RealVector _speed;

  /// array for temporary speed of sound
  RealVector _soundSpeed;

  /// array of temporary values
  RealMatrix _values;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;
  
  /// epsilon factor (must be >0)
  CFreal _dissipEps;
  
  /// flag telling not to recompute the artificial viscous coefficient mu_s
  CFuint _freeze_mu_s;
  
}; // end of class Euler2DCarbuncleFixSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Euler2DCarbuncleFixSourceTerm_hh
