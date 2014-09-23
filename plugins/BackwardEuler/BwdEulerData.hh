// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_BackwardEuler_BwdEulerData_hh
#define COOLFluiD_Numerics_BackwardEuler_BwdEulerData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a Data Object that is accessed by the different
  /// BackwardEulerCom 's that compose the BackwardEuler.
  /// @see BackwardEulerCom
  /// @author Andrea Lani
  /// @author Tiago Quintino
class BwdEulerData: public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  BwdEulerData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~BwdEulerData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre the pointer to LinearSystemSolver is not constant to
  ///      allow dynamic_casting
  void setLinearSystemSolver(Framework::MultiMethodHandle<Framework::LinearSystemSolver> lss)
  {
    m_lss = lss;
  }

  /// Get the linear system solver
  Framework::MultiMethodHandle<Framework::LinearSystemSolver> getLinearSystemSolver() const
  {
    cf_assert(m_lss.isNotNull());
    return m_lss;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BwdEuler";
  }

private: // data

  /// Linear system solver
  Framework::MultiMethodHandle<Framework::LinearSystemSolver> m_lss;

}; // end of class BwdEulerData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for BackwardEuler
typedef Framework::MethodCommand<BwdEulerData> BwdEulerCom;

/// Definition of a command provider for BackwardEuler
typedef Framework::MethodCommand<BwdEulerData>::PROVIDER BwdEulerComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_BackwardEuler_BwdEulerData_hh
