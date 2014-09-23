// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_ParJFSetupGMRESR_hh
#define COOLFluiD_Numerics_Petsc_ParJFSetupGMRESR_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header
#include "Petsc/BaseSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents
  *
  * @author Jiri Simonek
  *
  */
class ParJFSetupGMRESR : public BaseSetup {
public:

  /**
   * Constructor.
   */
  explicit ParJFSetupGMRESR(const std::string& name);
	
  /**
   * Destructor.
   */
  ~ParJFSetupGMRESR();
	
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
private: //helper functions
  
  /**
   * Set up the index mapping
   */
  void setIdxMapping();

  /**
   * Set up the matrix
   */
  void setMatrix(const CFuint localSize,
                 const CFuint globalSize);

  /**
   * Set up the vectors
   */
  void setVectors(const CFuint localSize,
                  const CFuint globalSize);
  /**
  * Set up the inner GMRES solver
  */
	//void setInnerSolver(const CFreal tol,
	//                   const CFuint maxIter,
	//                 const PCType preconditioner);

private:

  /// epsilon for computing numerical derivative
  CFreal _epsilon;

	/// Maximum number of GMRES iteration inside GMRESR solver
	//CFint _gmresRestart;

  /// Relative tolerance for GMRES inner solver
	//CFreal _gmresRelTol;

  /// Absolute tolerance for GMRES inner solver
	//CFreal _gmresAbsTol;

  /// Divergence tolerance for GMRES inner solver
	//CFreal _gmresDivTol;

  /// Preconditioner for GMRES inner solver
	//PCType _gmresPreconditioner;

  /// Maximum number of GMRESR iteration if it does not converge
  CFint _gmresrMaxIter;

  /// Tolerance for GMRESR solver
  CFreal _gmresrTol;

  /// Dimension of GMRESR solver subspace
  CFint _gmresrSubspaceDim;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_ParJFSetupGMRESR_hh
