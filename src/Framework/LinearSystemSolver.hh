// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LinearSystemSolver_hh
#define COOLFluiD_Framework_LinearSystemSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/Method.hh"
#include "Framework/MultiMethodHandle.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class LSSMatrix;
    class LSSVector;
    class BlockAccumulator;
    class LSSData;
    
//////////////////////////////////////////////////////////////////////////////

/// This class represents a LinearSystemSolver.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API LinearSystemSolver : public Method,
                           public Common::DynamicFunctionCaller<LinearSystemSolver> {

public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<LinearSystemSolver,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  LinearSystemSolver(const std::string& name);

  /// Default destructor
  virtual ~LinearSystemSolver();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Solve the linear system
  /// @post pushs and pops the Namespace to which this Method belongs
  void solveSys();

  /// Prints the Linear System to a file.
  virtual void printToFile(const std::string prefix, const std::string suffix) = 0;
  
  /// Create a block accumulator with chosen internal storage
  /// @return a newly created block accumulator
  /// @post the block has to be deleted outside
  virtual BlockAccumulator* createBlockAccumulator (const CFuint nbRows, 
						    const CFuint nbCols, 
						    const CFuint subBlockSize, 
						    CFreal* ptr = CFNULL) const = 0;
  
  /// Get the LSS idx mapping (from local to global LSS ids)
  /// Mutator that allows changes for clients.
  LSSIdxMapping& getLocalToGlobalMapping();
  
  /// Get the LSS idx mapping (from local to locally updatable LSS ids)
  /// -1 is the default value (meaning not updatable)
  /// Accessor that does not allow mutation for clients.
  LSSIdxMapping& getLocalToLocallyUpdatableMapping();
  
  /// Mask array that specifies which equations are solved by the current LSS
  Common::SafePtr<std::valarray<bool> > getMaskArray() {   return &m_maskArray;  }
  
  /// Array listing the IDs of all the equations to be solved
  /// by this LSS
  Common::SafePtr<std::vector<CFuint> > getEquationIDs() {  return &m_maskEquationIDs; }
  
  /// Gets the size of the system of equations to solve
  CFuint getNbSysEqs() const {   return m_nbSysEquations;  }

  /// Get the Preconditioner system matrix
  virtual Common::SafePtr<LSSMatrix> getPreconditionerMatrix() const
  {
    throw Common::NotImplementedException(FromHere(), "LinearSystemSolver::getPreconditionerMatrix()");
  }

  /// Get the LSS system matrix
  virtual Common::SafePtr<LSSMatrix> getMatrix() const = 0;

  /// Get the LSS right hand side vector
  virtual Common::SafePtr<LSSVector> getRhsVector() const = 0;

  /// Get the LSS solution vector
  virtual Common::SafePtr<LSSVector> getSolVector() const = 0;

  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<LinearSystemSolver>::run_dynamic_function(func);
  }

  /// Gets the Class name
  static std::string getClassName() {  return "LinearSystemSolver"; }

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

protected: // abstract interface implementations

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void unsetMethodImpl();

  /// Solve the linear system
  /// This is the abstract function that the concrete methods must implement.
  virtual void solveSysImpl() = 0;

  /// This mutator is meant to be used only by derived classes
  CFuint& getNbSysEquations() { return m_nbSysEquations;  }

protected: // method data
  
  /// the size of the system of equations to solve
  Common::SafePtr<LSSData> m_lssData;
  
  /// the size of the system of equations to solve
  CFuint m_nbSysEquations;

  /// mask array to specify which equations are solved by the current LSS
  std::valarray<bool> m_maskArray;

  /// array storing the IDs of the equations to solve with the current LSS
  std::vector<CFuint> m_maskEquationIDs;

}; // class LinearSystemSolver

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(LinearSystemSolver) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LinearSystemSolver_hh
