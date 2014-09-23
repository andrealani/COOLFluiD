// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_NewtonIteratorData_hh
#define COOLFluiD_Numerics_NewtonMethod_NewtonIteratorData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a Data Object that is accessed by the different
  /// NewtonMethodCom 's that compose the NewtonMethod.
  /// @see NewtonMethodCom
  /// @author Tiago Quintino
class NewtonIteratorData : public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  NewtonIteratorData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~NewtonIteratorData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NewtonIterator";
  }

  /// Gets the maximum number of steps to be performed
  /// @param iSys is the ID of the current subsystem of equations to solve
  CFuint getNbMaxSteps(CFuint iSys) const
  {
    cf_assert(iSys < m_maxSteps.size());
    return m_maxSteps[iSys];
  }

  /// Gets the maximum norm to be achieved
  CFreal getMaxNorm() const
  {
    return m_maxNorm;
  }

  /// Checks if convergence history should be printed.
  bool isPrintHistory() const
  {
    return m_printHistory;
  }

  /// Checks if linear system matrix, rhs and solution vectors should be saved to files
  bool isSaveSystemToFile() const
  {
    return m_saveSystemToFile;
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
  
private: // data

  /// flag to indicate that convergence has been achieved
  bool m_achieved;

  /// Linear system solver
  Framework::MultiMethodHandle<Framework::LinearSystemSolver> m_lss;

  /// number of newton iterations per global time step
  std::vector<CFuint> m_maxSteps;

  /// L2 norm of dU to reach per global time step
  CFreal m_maxNorm;

  /// flag to indicate printing of history in the newton iteration
  bool m_printHistory;

  /// flag to indicate saving files of system matrix, rhs and solution vectors at each iteration
  bool m_saveSystemToFile;

}; // end of class NewtonIteratorData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for NewtonMethod
typedef Framework::MethodCommand<NewtonIteratorData> NewtonIteratorCom;

/// Definition of a command provider for NewtonMethod
typedef Framework::MethodCommand<NewtonIteratorData>::PROVIDER NewtonIteratorComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_NewtonIteratorData_hh

