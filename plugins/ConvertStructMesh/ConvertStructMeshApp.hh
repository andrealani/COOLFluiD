// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshTools_ConvertStructMeshApp_hh
#define COOLFluiD_CFmeshTools_ConvertStructMeshApp_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class converts a structured mesh to unstructured Tecplot
 *
 * @author Andrea Lani
 *
 */
/**
 * This class defines the Module ConvertStructMeshApp
 */
class ConvertStructMeshAppModule : public Environment::ModuleRegister<ConvertStructMeshAppModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ConvertStructMeshApp";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements converters for structured and block structured meshes.";
  }

}; // end ConvertStructMeshApp

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshTools_ConvertStructMeshApp_hh
