// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_PetscOptions_hh
#define COOLFluiD_Numerics_Petsc_PetscOptions_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must be before any other header

#include "Common/StringOps.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a singleton storing mapping between
 * a string and the corresponding Petsc option
 *
 * @author Andrea Lani
 *
 */
class PetscOptions {
public:

  /**
   * Set all the options
   */
  static void setAllOptions();

  /**
   * Gets the PCType corresponding to a given string
   */
  static PCType getPCType(const std::string& name)
  {
    return _pcType.find(name)->second;
  }

  /**
   * Gets the KSPType corresponding to a given string
   */
  static KSPType getKSPType(const std::string& name)
  {
    return _kspType.find(name)->second;
  }

  /**
   * Gets the MatOrderingType corresponding to a given string
   */
  static MatOrderingType getMatOrderType(const std::string& name)
  {
    return _matOrderType.find(name)->second;
  }

private:
  /**
   * Default constructor without arguments.
   */
  PetscOptions();

  /**
   * Destructor.
   */
  ~PetscOptions();

  /**
   * Set all the PC types
   */
  static void setPCTypes();

  /**
   * Set all the KSP types
   */
  static void setKSPTypes();

  /**
   * Set all the mat odering types
   */
  static void setMatOrderTypes();

 private: // data

  /// mapping string - PCType
  static std::map<std::string, PCType, std::less<std::string> > _pcType;

  /// mapping string - KSPType
  static std::map<std::string, KSPType, std::less<std::string> > _kspType;

  /// mapping string - MatOrderingType
  static std::map<std::string,
		  MatOrderingType,
		  std::less<std::string> > _matOrderType;

}; // end of class PetscOptions

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_PetscOptions_hh
