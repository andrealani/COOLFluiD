// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_StdSeqSetup_hh
#define COOLFluiD_Numerics_Petsc_StdSeqSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Petsc/BaseSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class represents a command to be executed during the set up
  * of a the Petsc LSS method in sequential environments
  *
  * @author Andrea Lani
  * @author Tiago Quintino
  */
class StdSeqSetup : public BaseSetup {
public:

  /**
   * Constructor.
   */
  explicit StdSeqSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSeqSetup();

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

#endif // COOLFluiD_Numerics_Petsc_StdSeqSetup_hh

