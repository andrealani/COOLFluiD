// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_TwoLayerParSetup_hh
#define COOLFluiD_Numerics_Petsc_TwoLayerParSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Petsc/BaseSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a command to be executed during the set up
  * of a the Petsc LSS method in sequential environment, where two
  * sets of states are solved concurrently in the same jacobian matrix.
  *
  * @author Thomas Wuilbaut
  * @author Tiago Quintino
  *
  */
class TwoLayerParSetup : public BaseSetup {
public:

  /**
   * Constructor.
   */
  explicit TwoLayerParSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~TwoLayerParSetup();

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

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_TwoLayerParSetup_hh

