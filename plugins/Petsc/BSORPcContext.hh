// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_BSORPcContext_hh
#define COOLFluiD_Numerics_Petsc_BSORPcContext_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataStorage.hh"
#include "Petsc/PetscHeaders.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {
    class MFContext;
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a tuple of data to be passed to the BSOR
 * user defined preconditioner of Petsc
 *
 * @author Jiri Simonek 
 *
 */
class BSORPcContext {
public: // functions

  /// Constructor
  BSORPcContext() {}

  /// pointer to MFContext - pointer to preconditioner matrix is needed from MFContext
  MFContext* pMFC;

  /// Type of Gauss-Seidel relaxation
  //MatSORType relaxType;

}; // end of class BSORPcContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_BSORPcContext_hh
