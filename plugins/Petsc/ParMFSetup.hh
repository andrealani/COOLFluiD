// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_ParMFSetup_hh
#define COOLFluiD_Numerics_Petsc_ParMFSetup_hh

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
  * @author Jiri Simonek
  *
  */
class ParMFSetup : public BaseSetup {
public:

  /**
   * Constructor.
   */
  explicit ParMFSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~ParMFSetup();
	
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
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

private:

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_ParMFSetup_hh
