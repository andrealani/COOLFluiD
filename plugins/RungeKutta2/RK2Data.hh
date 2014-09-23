// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKutta2_RK2Data_hh
#define COOLFluiD_Numerics_RungeKutta2_RK2Data_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace RungeKutta2 {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// RungeKutta2Com 's that compose the RungeKutta2.
/// @see RungeKutta2Com
/// @author Tiago Quintino
class RK2Data : public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  RK2Data(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~RK2Data();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the flag that indicates the method is time accurate
  bool isTimeAccurate() const
  {
    return m_isTimeAccurate;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "RK2";
  }

private: // data

  /// flag to indicate if time accurate
  bool m_isTimeAccurate;

}; // end of class RK2Data

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for RungeKutta2
typedef Framework::MethodCommand<RK2Data> RK2Com;

/// Definition of a command provider for RungeKutta2
typedef Framework::MethodCommand<RK2Data>::PROVIDER RK2ComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta2

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKutta2_RK2Data_hh
