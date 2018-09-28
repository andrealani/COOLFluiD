// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_FwdEulerData_hh
#define COOLFluiD_Numerics_ForwardEuler_FwdEulerData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "MathTools/RealVector.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"

#include "ForwardEuler/ForwardEulerAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a Data Object that is accessed by the different
  /// ForwardEulerCom 's that compose the ForwardEuler.
  /// @see ForwardEulerCom
  /// @author Tiago Quintino
class  ForwardEuler_API FwdEulerData : public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  FwdEulerData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~FwdEulerData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FwdEuler";
  }

  /// Sets the current norm
  void setNorm(CFreal p)
  {
    m_norm[0] = p;
  }

  /// Sets the current norm
  void setNorm(RealVector p)
  {
    cf_assert(m_norm.size() == p.size());
    m_norm = p;
  }

  /// Sets the size of m_norm vector
  void setNormSize(CFuint p)
  {
    m_norm.resize(p);
  }

  /// Gets the current norm
  RealVector getNorm() const
  {
    return m_norm;
  }

  /// Gets the variable in which to compute the norm
  CFuint getVarID() const
  {
    return m_varID;
  }


  /// Checks if convergence history should be printed.
  bool isPrintHistory() const
  {
    return m_printHistory;
  }

  /// Gets the flag that indicates the method is time accurate
  bool isTimeAccurate() const
  {
    return m_isTimeAccurate;
  }
  
  /// Gets the flag that indicates if global time step is used
  bool isGlobalTimeStep() const
  {
    return m_isGlobalTimeStep;
  }
  
  /// Gets the flag that indicates we are at the last iteration
  bool isAchieved() const
  {
    return m_achieved;
  }

  /// Sets the flag that indicates we are at the last iteration
  void setAchieved(bool achieved)
  {
    m_achieved = achieved;
  }

private: // data

  /// flag to indicate that convergence has been achieved
  bool m_achieved;

  /// L2 norm of dU
  RealVector m_norm;

  /// Variable for the L2 norm computation
  CFuint m_varID ;

  /// flag to indicate printing of history of iterations
  bool m_printHistory;

  /// flag to indicate if time accurate
  bool m_isTimeAccurate;

  /// flag to indicate if global time stepping
  bool m_isGlobalTimeStep;
  
}; // end of class FwdEulerData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for ForwardEuler
typedef Framework::MethodCommand<FwdEulerData> FwdEulerCom;

/// Definition of a command provider for ForwardEuler
typedef Framework::MethodCommand<FwdEulerData>::PROVIDER FwdEulerComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_FwdEulerData_hh
