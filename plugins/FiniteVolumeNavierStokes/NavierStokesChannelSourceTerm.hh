#ifndef COOLFluiD_Numerics_FiniteVolume_NavierStokesChannelSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_NavierStokesChannelSourceTerm_hh

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
 * This class represents a source term for a channel flow.
 * It simulates a pressure gradient in the streamwise direction, while
 * the real pressure in the channel remains constant. This is useful
 * to use periodic boundary conditions connecting inlet and outlet.
 * 
 * Theory: dp/dx = - (u*)^2 rho / delta
 *  where u* = friction velocity
 *        rho = density
 *        delta = half channel width
 * 
 * The term -(dp/dx) is added to the RHS as source
 * 
 * @author Willem Deconinck
 *
 */
template <class EULERVAR, class NSVAR>
class NavierStokesChannelSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  NavierStokesChannelSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesChannelSourceTerm();

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
  
  
  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "NavierStokesChannelSourceTerm";
  }
  
protected: // data

  /// Volumes data handle
  Framework::DataHandle<CFreal> m_volumes;
  
  /// corresponding variable set
  Common::SafePtr<EULERVAR> m_varSet;

  /// corresponding diffusive variable set
  Common::SafePtr<NSVAR> m_diffVarSet;
  
  /// Euler physical data
  RealVector m_physicalData;

  /// ID for the streamwise momentum equation where the source term will be applied to
  CFuint m_streamwiseMomentumEquationID;

  /// friction velocity
  CFreal m_frictionVelocity;
  
  /// half channel width
  CFreal m_halfChannelWidth;

  /// Density in source term
  CFreal m_rho;

  /// flag if density is variable or constant
  bool m_useVariableDensity;
  
  /// channel length (for informative calculation of pressure drop over the channel, optional)
  CFreal m_channelLength;
  
}; // end of class NavierStokesChannelSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesChannelSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokesChannelSourceTerm_hh
