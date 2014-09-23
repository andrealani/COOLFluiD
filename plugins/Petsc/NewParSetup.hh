// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_NewParSetup_hh
#define COOLFluiD_Numerics_Petsc_NewParSetup_hh

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
  * @author Andrea Lani
  *
  */
class NewParSetup : public BaseSetup {
public:

  /**
   * Constructor.
   */
  explicit NewParSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~NewParSetup();

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


}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_NewParSetup_hh
