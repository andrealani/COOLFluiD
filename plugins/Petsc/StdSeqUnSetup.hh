// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_StdSeqUnSetup_hh
#define COOLFluiD_Numerics_Petsc_StdSeqUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Petsc/PetscLSSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Andrea Lani
 */
class StdSeqUnSetup : public PetscLSSCom {
public:

  /**
   * Constructor.
   */
  explicit StdSeqUnSetup(const std::string& name) : PetscLSSCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdSeqUnSetup()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

}; // class UnSeqSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_StdSeqUnSetup_hh

