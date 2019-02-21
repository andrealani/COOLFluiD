// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConvergenceMethodData_hh
#define COOLFluiD_Framework_ConvergenceMethodData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FilterState.hh"
#include "Framework/FilterRHS.hh"
#include "Framework/MethodData.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/CFL.hh"
#include "Framework/SpaceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
   class VarRegistry;
 
//////////////////////////////////////////////////////////////////////////////

/// Base class for the Data of the ConvergenceMethod's.
/// Provides the computation of the norm.
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API ConvergenceMethodData: public Framework::MethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ConvergenceMethodData(Common::SafePtr<Method> owner);

  /// Destructor.
  virtual ~ConvergenceMethodData();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Get the Norm computer
  Common::SafePtr<Framework::ComputeNorm> getNormComputer() const { return m_computeNorm.getPtr(); }

  /// Get the CFL
  Common::SafePtr<Framework::CFL> getCFL() { return &m_CFL; }

  /// Get an object that is able to filter the state values
  /// For example, if pressure is below zero, set to zero.
  /// @returns a FilterState to be used by the concrete methods for filtering
  Common::SafePtr<Framework::FilterState> getFilterState() const {  return m_filterState.getPtr(); }

  /// Get an object that is able to filter the RHS values
  /// @returns a FilterRHS to be used by the concrete methods for filtering
  Common::SafePtr<Framework::FilterRHS> getFilterRHS() const {  return m_filterRHS.getPtr(); }
  
  /// Updates the residual and places it in the SubSystemStatus
  void updateResidual();

  /// Get the solving rate
  CFuint getSolvingRate() const
  {
    return m_solvingRate;
  }
  
  /// Reset the flag telling if to compute the jacobian
  void setDoComputeJacobFlag(bool doComputeJacobian)
  {
    m_doComputeJacob = doComputeJacobian;
  }
  
  /// Get the flag telling if to compute the jacobian
  bool getDoComputeJacobFlag() const
  {
    return m_doComputeJacob;
  }

  /// Get the flag telling if to update the solution
  bool getDoUpdateSolution() const
  {
    return m_doUpdateSolution;
  }
  
  /// Get the flag telling if to freeze the jacobian
  bool freezeJacobian() const
  {
    return m_freezeJacobian;
  }
  
  /// Flag telling whether only preprocessing the solution once
  bool onlyPreprocessSolution() const
  {
    return m_onlyPreprocessSolution;
  }

  /// Get the number of systems to solve in the convergence method
  CFint getNbLSSToSolveAtOnce() const
  {
    return m_nbLSSToSolveAtOnce;
  }
  
  /// Get the flag telling if to update the solution
  ConvergenceStatus& getConvergenceStatus()
  {
    return m_cstatus;
  }
  
  /// Compute and store the space residual norm
  void computeSpaceResidualNorm();
  
  /// Get the space residual
  RealVector getSpaceResidual() const ;
  
  /// Set  the space residual
  void setSpaceResidual(const RealVector& spaceResidual);
   
  /// Sets the SpaceMethod for this SpaceMethod to use
  /// @pre the pointer to SpaceMethod is not constant to
  ///      allow dynamic_casting
  void setSpaceMethod(Framework::MultiMethodHandle<Framework::SpaceMethod> sm) 
  {
    m_sm = sm;
  }

  /// Get the space method
  Framework::MultiMethodHandle<Framework::SpaceMethod> getSpaceMethod() const
  {
    cf_assert(m_sm.isNotNull());
    return m_sm;
  }
  
  /// Set the factory registry
  virtual void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr);
  
private: // data

  /// Current convergence status
  Framework::ConvergenceStatus m_cstatus;
  
  /// Space Method
  Framework::MultiMethodHandle<Framework::SpaceMethod> m_sm;
  
  /// CFL object
  Framework::CFL  m_CFL;
    
  /// Name of the norm
  std::string m_normStr;

  /// Functor that computes the requested norm
  Common::SelfRegistPtr<Framework::ComputeNorm>  m_computeNorm;

  /// Name of the filter state object
  std::string m_filterStateStr;

  /// Name of the filter RHS object
  std::string m_filterRHSStr;
  
  /// Filter state to apply some fixes
  Common::SelfRegistPtr<Framework::FilterState>  m_filterState;
  
  /// Filter RHS to apply some fixes
  Common::SelfRegistPtr<Framework::FilterRHS>  m_filterRHS;

  /// flag telling the solving rate
  CFuint m_solvingRate;
  
  /// flag to tell to compute the jacobian
  bool m_doComputeJacob;
  
  /// flag to tell to update the solution
  bool m_doUpdateSolution;
  
  /// flag to tell to freeze the jacobian during the iterative process
  bool m_freezeJacobian;

  /// flag to tell to preprocess the solution once
  bool m_onlyPreprocessSolution;
    
  /// number of linear system solvers to solve at once
  CFint m_nbLSSToSolveAtOnce;
  
  Common::SafePtr<VarRegistry> ssys_var_regist;
  
}; // end of class ConvergenceMethodData

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvergenceMethodData_hh
