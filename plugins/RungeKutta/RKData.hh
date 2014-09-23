// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKutta_RKData_hh
#define COOLFluiD_Numerics_RungeKutta_RKData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// RungeKuttaCom 's that compose the RungeKutta.
/// @see RungeKuttaCom
/// @author Thomas Wuilbaut
class RKData : public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  RKData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~RKData();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the flag that indicates if we are at the first step
  bool getIsFirstStep() const
  {
    return m_isFirst;
  }

  /// Sets the flag that indicates if we are at the first step
  void setIsFirstStep(const bool first)
  {
    m_isFirst = first;
  }

  /// Gets the flag that indicates if we are at the last step
  bool getIsLastStep() const
  {
    return m_isLast;
  }

  /// Sets the flag that indicates if we are at the last step
  void setIsLastStep(const bool last)
  {
    m_isLast = last;
  }

  /// Gets the alpha coefficient
  CFreal getAlpha(const CFuint i) const
  {
    return m_alpha[i];
  }

  /// Gets the beta coefficient
  CFreal getBeta(const CFuint i) const
  {
    return m_beta[i];
  }

  /// Gets the gamma coefficient
  CFreal getGamma(const CFuint i) const
  {
    return m_gamma[i];
  }

  /// Gets the order of the method
  CFuint getOrder() const
  {
    return m_order;
  }

  /// Gets the current step
  CFuint getCurrentStep() const
  {
    return m_step;
  }

  /// Sets the current step
  void setCurrentStep(const CFuint step)
  {
    m_step = step;
  }

  /// Gets the flag that indicates the method is time accurate
  bool isTimeAccurate() const
  {
    return m_isTimeAccurate;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "RKLS";
  }

  /// @return m_timeIterN
  CFreal getTimeIterN()
  {
    return m_timeIterN;
  }

  /// sets m_timeIterN
  void setTimeIterN(const CFreal timeIterN)
  {
    m_timeIterN = timeIterN;
  }

private:

  /// flag to check first iteration
  bool m_isFirst;

  /// flag to check last iteration
  bool m_isLast;

  /// Iterator used for counting the RK method steps
  CFuint m_step ;

  /// Order of the RK method
  CFuint m_order;

  /// vector of the coeficients of the RK method
  std::vector<CFreal> m_alpha;

  /// vector of the coeficients of the RK method
  std::vector<CFreal> m_beta;

  /// vector of the time coeficients of the RK method
  std::vector<CFreal> m_gamma;

  /// flag to indicate if time accurate
  bool m_isTimeAccurate;

  /// time at beginning of iteration
  CFreal m_timeIterN;

}; // end of class RKData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for RungeKutta
typedef Framework::MethodCommand<RKData> RKCom;

/// Definition of a command provider for RungeKutta
typedef Framework::MethodCommand<RKData>::PROVIDER RKComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKutta_RKData_hh
