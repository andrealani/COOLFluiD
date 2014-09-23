#ifndef COOLFluiD_Numerics_FiniteVolume_NavierStokes2DAxiSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_NavierStokes2DAxiSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
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
template <class EULERVAR, class NSVAR>
class NavierStokes2DAxiSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  NavierStokes2DAxiSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokes2DAxiSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
    
    _sockets.template createSocketSink<RealVector>("nstates");
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
    
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<EULERVAR> _varSet;

  /// corresponding diffusive variable set
  Common::SafePtr<NSVAR> _diffVarSet;

  /// vector to store temporary result
  RealVector _temp;
  
  /// Euler physical data
  RealVector _physicalData;
  
  /// pressure jacobian
  RealVector _pressureJacob;

  /// array of temporary values
  RealMatrix _values;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;

  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;
  
  /// vector of IDs for u and v components
  std::vector<CFuint> _uvID;
  
  /// use the gradient computed with the least square reconstruction
  bool _useGradientLS;
  
}; // end of class NavierStokes2DAxiSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DAxiSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokes2DAxiSourceTerm_hh
