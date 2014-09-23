#ifndef COOLFluiD_AeroCoef_AeroCoefStopCondition_hh
#define COOLFluiD_AeroCoef_AeroCoefStopCondition_hh

//////////////////////////////////////////////////////////////////////////////

#include <deque>

#include "Framework/StopCondition.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a StopConditon for the convergence process
 * based on the convergence of CL, CD and CM.
 *
 * @author Tiago Quintino
 */
class AeroCoefStopCondition : public Framework::StopCondition {

private: // class

struct StoreCoefs
{
  /// storage of CL
  CFreal CL;
  /// storage of CD
  CFreal CD;
  /// storage of CM
  CFreal CM;
};

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  AeroCoefStopCondition(const std::string& name);

  /// Default destructor
  virtual ~AeroCoefStopCondition();

  /// Configure this object with user defined parameters
  /// @param args arguments for the configuration
  virtual void configure ( Config::ConfigArgs& args );

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

  /// number of iterations to use as information to compute the derivatives
  CFuint m_nb_iter;

 /// maximum number of iterations 
  CFuint m_max_iter;

  /// converge threshold in CL
  CFreal m_conv_cl;

  /// converge threshold in CD
  CFreal m_conv_cd;

  /// converge threshold in CM
  CFreal m_conv_cm;

  /// check convergence of CL
  bool m_check_cl;

  /// check convergence of CD
  bool m_check_cd;

  /// check convergence of CM
  bool m_check_cm;

  /// storage of the values in the last iterations
  std::deque<StoreCoefs> m_past_coeffs;

}; // end of class AeroCoefStopCondition

//////////////////////////////////////////////////////////////////////////////

  } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_AeroCoef_AeroCoefStopCondition_hh
