// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_ILUPreconditioner_hh
#define COOLFluiD_Numerics_Petsc_ILUPreconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "Petsc/ShellPreconditioner.hh"
#include "Petsc/ILUPcContext.hh"

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
 * @author Jiri Simonek
 *
 */

class ILUPreconditioner : public ShellPreconditioner {
public:

  /**
   * Constructor
   */
  ILUPreconditioner(const std::string& name);

  /**
   * Default destructor
   */
  ~ILUPreconditioner();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

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

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// ILU context 
  ILUPcContext _pcc;

  /// Petsc Index Set
  IS _indexSet;

  /// Petsc Index Set Array
  CFint* _indexSetArray;

  /// Petsc Matrix Factorization Info
  MatFactorInfo _info;

}; // end of class ILUPreconditioner

//////////////////////////////////////////////////////////////////////////////

  } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_ILUPreconditioner_hh
