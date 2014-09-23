// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_hh
#define COOLFluiD_Numerics_Trilinos_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement the interface to Trilinos Linear System Solver
  namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

      /// preconditioner type
      typedef int PCType;

      /// domain decomposition sub solve type
      typedef int PCSubSolveType;

      /// krylov solver type
      typedef int KSPType;

      /// residual type
      typedef int ResidualType;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Trilinos
 */
class TrilinosModule : public Environment::ModuleRegister<TrilinosModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Trilinos";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Trilinos linear system solver.";
  }

}; // end TrilinosModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_Trilinos_hh

