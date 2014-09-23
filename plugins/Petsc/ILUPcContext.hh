// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_ILUPcContext_hh
#define COOLFluiD_Numerics_Petsc_ILUPcContext_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataStorage.hh"
#include "Common/ConnectivityTable.hh"
#include "PetscMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {
    class JFContext;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a tuple of data to be passed to the (Block)ILU preconditioner for Jacobian-Free
 * user defined preconditioner of Petsc
 *
 * @author Jiri Simonek 
 *
 */
class ILUPcContext {
public: // functions

  /// Constructor
  ILUPcContext() {}

  /// pointer to JFContext
  JFContext* pJFC;

  /// preconditioner matrix
  PetscMatrix* precondMat;

}; // end of class ILUPcContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_ILUPcContext_hh
