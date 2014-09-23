// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_NullPreconditioner_hh
#define COOLFluiD_Numerics_Petsc_NullPreconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "Petsc/ShellPreconditioner.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {
    class MatrixInverter;
  }

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a shell preconditioner object
 *
 * @author Andrea Lani
 *
 */

class NullPreconditioner : public ShellPreconditioner {
public:

  /**
   * Constructor
   */
  NullPreconditioner(const std::string& name);

  /**
   * Default destructor
   */
  ~NullPreconditioner();

  /**
   * Set the preconditioner
   */
  virtual void setPreconditioner();

  /**
   * Compute before solving the system
   */
  virtual void computeBeforeSolving();
  
  /**
   * Compute after solving the system
   */
  virtual void computeAfterSolving();
  
}; // end of class NullPreconditioner

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_NullPreconditioner_hh
