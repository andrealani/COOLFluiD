#ifndef COOLFluiD_SubSystemCoupler_HeatFluxStopCondition_hh
#define COOLFluiD_SubSystemCoupler_HeatFluxStopCondition_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StopCondition.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

      namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a StopConditon for the convergence process
 * based on the convergence of the L2norm concerning the difference of the fluxes
 * on solid and fluid common interface.
 *
 * @author Leonardo Nettis
 */
class HeatFluxStopCondition : public  Framework::StopCondition {

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  HeatFluxStopCondition(const std::string& name);

  /// Default destructor
  virtual ~HeatFluxStopCondition();

  /// Configure this object with user defined parameters
  /// @param args arguments for the configuration
  virtual void configure (const Config::ConfigArgs& args);

  /**
   * returns true if we have to evaluate the stopcondition
   * needs to combine values from different CPU domains
   */
  virtual bool IsGlobal () const;

  /**
   * Take the combined value from all CPU's and decide if the SubSystem
   * should stop
   */
  virtual bool isAchieved (const Framework::ConvergenceStatus& status) ;

private: // data

  /// converge threshold
  CFreal m_conv_norm;

  /// norm value on first iteration
  CFreal m_first_norm;

}; // end of class HeatFluxStopCondition

//////////////////////////////////////////////////////////////////////////////

      } // end of namespace SubSystemCoupler

  } // end of namespace Numerics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SubSystemCoupler_HeatFluxStopCondition_hh


