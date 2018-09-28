// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SpaceMethodData_hh
#define COOLFluiD_Framework_SpaceMethodData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/MethodData.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/NumericalJacobian.hh"
#include "Framework/LinearSystemSolver.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

        class CFL;
        class ConvergenceMethod;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a data object that is accessed by the different
/// SpaceMethodCom 's that compose the SpaceMethod.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API SpaceMethodData : public Framework::MethodData {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  SpaceMethodData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  virtual ~SpaceMethodData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Get the solution variables set
  Common::SafePtr<Framework::ConvectiveVarSet> getSolutionVar() const
  {
    cf_assert(_solutionVar.isNotNull());
    return _solutionVar.getPtr();
  }

  /// Get the update variables set
  Common::SafePtr<Framework::ConvectiveVarSet> getUpdateVar() const
  {
    cf_assert(_updateVar.isNotNull());
    return _updateVar.getPtr();
  }

  /// Get the diffusive variables set
  Common::SafePtr<Framework::DiffusiveVarSet> getDiffusiveVar() const
  {
    cf_assert(_diffusiveVar.isNotNull());
    return _diffusiveVar.getPtr();
  }

  /// Get the update variables set name
  std::string getUpdateVarStr() const
  {
    return _updateVarStr;
  }

  ///  Get the solution variables set name
  std::string getSolutionVarStr() const
  {
    return _solutionVarStr;
  }

  /// Get the diffusive variables set name
  std::string getDiffusiveVarStr() const
  {
    return _diffusiveVarStr;
  }
  
  class PreconditionerData {
  public:
    CFuint currentStateID;
    bool useBiggerStateIDs;
    bool useAllStateIDs;
    RealVector result;
  };
  
  /// set is restart flag
  void setIsRestart(bool isRestart) { m_isRestart = isRestart;}	

  /// get is restart flag
  bool isRestart() const { return m_isRestart;}
 
  /**
   * Get the numerical jacobian calculator
   */
  Framework::NumericalJacobian& getNumericalJacobian()
  {
    return _numericalJacobian;
  }
  
  /// preconditioner data
  Common::SafePtr<SpaceMethodData::PreconditionerData> getPreconditionerData()
  {
    return &_preconditionerData;
  }
      
  /// Get the ConvergenceMethod
  virtual Framework::MultiMethodHandle<Framework::ConvergenceMethod> getConvergenceMethod() const = 0;
  
  /**
   * Sets the LinearSystemSolver for this SpaceMethod to use
   * @pre the pointer to LinearSystemSolver is not constant to
   *      allow dynamic_casting
   */
  void setLinearSystemSolver(Framework::MultiMethodHandle<Framework::LinearSystemSolver> lss)
  {
    _lss = lss;
  }

  /**
   * Get the linear system solver
   */
  Framework::MultiMethodHandle<Framework::LinearSystemSolver> getLinearSystemSolver() const
  {
    cf_assert(_lss.isNotNull());
    return _lss;
  }
  
  /// Get the CFL from the convergence method
  Common::SafePtr<Framework::CFL> getCFL();

  /// Set the flag telling if to preprocess solution only once
  void setOnlyPreprocessSolution(bool onlyPreprocessSolution)
  {
    _onlyPreprocessSolution = onlyPreprocessSolution;
  }
  
  /// Get the flag telling if to preprocess solution only once
  bool doOnlyPreprocessSolution() const
  {
    return _onlyPreprocessSolution;
  }
  
  /// Set the flag telling if the jacobian has to be computed
  void setComputeJacobianFlag(bool computeJacobian)
  {
    _computeJacobian = computeJacobian;
  }

  /// Get the flag telling if the jacobian has to be computed
  bool doComputeJacobian() const
  {
    return _computeJacobian;
  }

  /// Check if the system matrix is frozen. If not, it should be recomputated on each iteration.
  bool isSysMatrixFrozen() const
  {
    return _sysMatFrozen;
  }

  /// Set the freezed state of the system matrix
  void setSysMatrixFrozen(bool frozen)
  {
    _sysMatFrozen = frozen;
  }

  /// Check if the system matrix will be frozen every iteration
  bool isSysMatrixFreezedEveryIteration() const
  {
    return _freezeSysMatEverIter;
  }
  
  /// get flag telling if the flux is being perturbed
  bool isPerturb() const {return _isPerturb;}
  
  /// set flag telling if the flux is being perturbed
  void setIsPerturb(bool isPerturb) {_isPerturb = isPerturb;}
  
  /// get flag telling the variables which is perturbed
  CFuint iPerturbVar() const {return _iPerturbVar;}
  
  /// set flag telling the variables which is perturbed
  void setIPerturbVar(CFuint iPerturbVar) {_iPerturbVar = iPerturbVar;}
  
  /// Flag telling to fill the preconditioner matrix
  bool fillPreconditionerMatrix() const
  {
    return _fillPreconditionerMatrix;
  }
  
  /// Flag telling to fill the preconditioner matrix
  void setFillPreconditionerMatrix(bool flag)
  {
    _fillPreconditionerMatrix = flag;
  }
  
  /**
   * Get the linear system matrix
   */
  Common::SafePtr<LSSMatrix> getLSSMatrix(CFuint iSys)
  {
    return (!_fillPreconditionerMatrix) ? _lss[iSys]->getMatrix() : _lss[iSys]->getPreconditionerMatrix();
  }
  
protected:

  /// Solution variables set
  Common::SelfRegistPtr<Framework::ConvectiveVarSet> _solutionVar;

  /// Update variables set
  Common::SelfRegistPtr<Framework::ConvectiveVarSet> _updateVar;

  /// Diffusive variables set
  Common::SelfRegistPtr<Framework::DiffusiveVarSet> _diffusiveVar;

  /// Linear system solver
  Framework::MultiMethodHandle<Framework::LinearSystemSolver> _lss;
  
  // the numerical jacobian computer
  Framework::NumericalJacobian _numericalJacobian;
  
  /// preconditioner data object
  PreconditionerData _preconditionerData;
  
  /// flag telling if the flux is being perturbed
  bool _isPerturb;
  
  /// flag telling if this is a restart
  bool m_isRestart;  
 
  /// index telling the variable being perturbed
  CFuint _iPerturbVar;
  
  /// Flag to fill the preconditioner matrix
  bool _fillPreconditionerMatrix;

  /// flag telling if to preprocess solution only once
  bool _onlyPreprocessSolution;
  
  /// flag telling if the jacobian has to be computed
  bool _computeJacobian;
  
  /// Freezed System matrix
  bool _sysMatFrozen;

  /// Strategy to freeze system matrix
  bool _freezeSysMatEverIter;
    
  /// string for configuration of the update variables
  std::string _updateVarStr;

  /// solution variable set name
  std::string _solutionVarStr;

  /// diffusive variable set name
  std::string _diffusiveVarStr;

}; // end of class SpaceMethodData

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SpaceMethodData_hh
