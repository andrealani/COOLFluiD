// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_AnalyticalModel_hh
#define COOLFluiD_Numerics_AnalyticalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement an DomainModel which provides
  /// an analytical description of the computational domain
  namespace AnalyticalModel {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module AnalyticalModel
 *
 * @author Tiago Quintino
 */
class AnalyticalModelModule : public Environment::ModuleRegister<AnalyticalModelModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "AnalyticalModel";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an analytical description of the computational domain";
  }

}; // end AnalyticalModelModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalModel

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_AnalyticalModel_hh

