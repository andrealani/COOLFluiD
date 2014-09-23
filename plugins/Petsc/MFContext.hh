// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_MFContext_hh
#define COOLFluiD_Numerics_Petsc_MFContext_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NumericalCommand.hh"
#include "Framework/DataStorage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {
    class PetscLSSData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a tuple of data to be passed to the matrix free solver of
 * Petsc
 *
 * @author Jiri Simonek
 *
 */
class MFContext {
public: // functions

  /// Constructor
  MFContext() {}

  /// pointer to the Petsc method data
  Common::SafePtr<PetscLSSData> petscData;

  /// Linear system matrix
  Common::SafePtr<PetscMatrix> mat;

}; // end of class MFContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_MFContext_hh
