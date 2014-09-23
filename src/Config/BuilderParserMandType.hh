// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_BuilderParserMandType_h
#define COOLFluiD_Config_BuilderParserMandType_h

/////////////////////////////////////////////////////////////////////////////

#include "Config/Config.hh"

namespace COOLFluiD
{
 namespace Config
 {

/////////////////////////////////////////////////////////////////////////////

  /// @brief Defines available mandatoriness policies.
  
  /// @author Quentin Gasper
  
  enum BuilderParserMandType
  {
   /// @brief Undefined policy
   MAND_UNDEFINED,
   
   /// @brief The presence of an element is forbidden
   MAND_FORBIDDEN,
   
   /// @brief The presence of an element is optional.
   
   /// The element may be not present or present with an empty value.
   MAND_OPTIONAL,
   
   /// @brief The presence of an element is mandatory.
   
   /// The element must be present with a non-empty value.
   MAND_MANDATORY
  }; // enum BuilderParserMandType

/////////////////////////////////////////////////////////////////////////////

 } // namespace Config
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_BuilderParserMandType_h