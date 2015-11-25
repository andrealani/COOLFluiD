// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_PetscLSSData_hh
#define COOLFluiD_Numerics_Petsc_PetscLSSData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Framework/LSSData.hh"

#include "Petsc/PetscVector.hh"
#include "Petsc/PetscMatrix.hh"
#include "Petsc/PetscOptions.hh"
#include "Petsc/JFContext.hh"
#include "Petsc/MFContext.hh"
#include "Petsc/ShellPreconditioner.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a data object that is accessed by the different
* PetscLSSCom 's that compose the PetscLSS.
 *
 * @todo there is missing documentation in this class.
 *
 * @author Andrea Lani
 * @author Jiri Simonek
 *
 */
class PetscLSSData : public Framework::LSSData {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  PetscLSSData(Common::SafePtr<std::valarray<bool> > maskArray,
	       CFuint& nbSysEquations,
	       Common::SafePtr<Framework::Method> owner);

  /**
   * Destructor.
   */
  ~PetscLSSData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();
  
  /**
   * Gets the shell preconditioner
   */
  Common::SafePtr<ShellPreconditioner> getShellPreconditioner()
  {
    cf_assert(_shellPreco.isNotNull());
    return _shellPreco.getPtr();
  }
  
  /**
   * Gets the solution vector
   */
  PetscVector& getSolVector()
  {
    return _xVec;
  }

  /**
   * Gets the rhs vector
   */
  PetscVector& getRhsVector()
  {
    return _bVec;
  }

  /**
   * Gets the matrix
   */
  PetscMatrix& getMatrix()
  {
    return _aMat;
  }

  /**
   * Gets the preconditioner matrix
   */
  PetscMatrix& getPreconditionerMatrix()
  {
    return _aPrecoMat;
  }

  /**
   * Gets the Jacobian-Free matrix
   */
  PetscMatrix& getJFMatrix()
  {
    return _jfMat;
  }

  /**
   * Context for the matrix free solver
   */
  JFContext* getJFContext()
  {
    return &_jfContext;
  }
  
  /**
   * Context for the matrix free solver
   */
  MFContext* getMFContext()
  {
    return &_mfContext;
  }
  
  /**
   * Gets the Preconditioner
   */
  PC& getPreconditioner()
  {
    return _pc;
  }

  /**
   * Gets the Krylov subspace method
   */
  KSP& getKSP()
  {
    return _ksp;
  }

  /**
   * Preconditioner type
   */
  PCType getPCType() const
  {
    return PetscOptions::getPCType(_pcTypeStr);
  }

  /**
   * Krylov solver type
   */
  KSPType getKSPType() const
  {
    return PetscOptions::getKSPType(_kspTypeStr);
  }

  /**
   * Matrix ordering type
   */
  MatOrderingType getMatOrderType() const
  {
    return PetscOptions::getMatOrderType(_matOrderTypeStr);
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "PetscLSS";
  }

  /**
   * Gets the relative tolerance
   */
  CFreal getRelativeTol() const
  {
    return _rTol;
  }

  /**
   * Gets the absolute tolerance
   */
  CFreal getAbsoluteTol() const
  {
    return _aTol;
  }

  /**
   * Gets the Divergence tolerance
   */
  CFreal getDivergenceTol() const
  {
    return _dTol;
  }

  /**
   * Show the Matrix Structure before solving
   */
  bool isShowMatrixStructure() const
  {
    return _showMatrixStructure;
  }

  /**
   * KSP convergence show rate
   */
  CFuint getKSPConvergenceShowRate() const
  {
    return _kspConvergenceShowRate;
  }
  
  /**
   * Gets the number of krylov subspaces
   */
  CFuint getNbKSP() const
  {
    return _nbKsp;
  }

  /**
   * Gets the number of levels of fill for the ILU preconditioner
   */
  CFuint getILULevels() const
  {
    //    return (_pcTypeStr == "PCILU"? _ilulevels:0);
    return _ilulevels;
  }

  /**
   * Get the info if user wants to use different matrix for preconditioner
   */
  bool getDifferentPreconditionerMatrix()
  {
    return _differentPreconditionerMatrix;
  }
  
  /**
   * Tell if AIJ structure must be used
   */
  bool useAIJ() {return _useAIJ;}
  
private:

  /// Shell preconditioner
  Common::SelfRegistPtr<ShellPreconditioner> _shellPreco;

  /// solution vector
  PetscVector _xVec;

  /// rhs vector
  PetscVector _bVec;
  
  /// matrix
  PetscMatrix _aMat;
  
  /// preconditioner matrix
  PetscMatrix _aPrecoMat;

  /// Jacobian-free matrix
  PetscMatrix _jfMat;

  /// Preconditioner
  PC _pc;

  /// Krylov subspace method
  KSP _ksp;

  /// jacobian free context
  JFContext _jfContext;

  /// jacobian free context2
  MFContext _mfContext;

  /// Number of Krylov spaces
  CFuint _nbKsp;

  /// Number of levels of fill for the ILU preconditioner
  CFuint _ilulevels;
  
  /// KSP convergence show rate
  CFuint _kspConvergenceShowRate;
  
  /// Preconditioner type
  std::string _pcTypeStr;

  /// Krylov solver type
  std::string _kspTypeStr;

  /// Mat ordering type
  std::string _matOrderTypeStr;

  /// name of the shell preconditioner
  std::string _shellPrecoStr;

  /// relative tolerance for iterative solver
  CFreal _rTol;

  /// absolute tolerance for iterative solver
  CFreal _aTol;

  /// divergence tolerance for iterative solver
  CFreal _dTol;

  ///flag if to show the matrix structure in X window
  bool _showMatrixStructure;

  /// Enable/disable usage of different matrix for preconditioner
  bool _differentPreconditionerMatrix;

  /// Use the AIJ structure instead of BAIJ
  bool _useAIJ;
    
}; // end of class PetscLSSData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for Petsc
typedef Framework::MethodCommand<PetscLSSData> PetscLSSCom;

/// Definition of a command provider for Petsc
typedef Framework::MethodCommand<PetscLSSData>::PROVIDER PetscLSSComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_PetscLSSData_hh
