// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_BuilderParserFrameInfo_h
#define COOLFluiD_Config_BuilderParserFrameInfo_h

//////////////////////////////////////////////////////////////////////////////

#include <string>
#include <map>

#include "Config/Config.hh"
#include "Common/xmlParser.h"
#include "Common/COOLFluiD.hh"

namespace COOLFluiD
{
 namespace Config
 {

//////////////////////////////////////////////////////////////////////////////

  /// @brief Holds a frame imformation
  
  struct Config_API BuilderParserFrameInfo
  {
   public:
    
    /// @brief Frame type    
    unsigned int frameType;
    
    /// @brief Frame attributes
    
    /// Key is the attribute names. Value is the attribute values
    std::map<std::string, std::string> frameAttributes;
    
    /// @brief Frame data
    XMLNode frameData;
    
    /// @brief Constructor
    BuilderParserFrameInfo();
    
    /// @brief Sets frame type.
    
    /// Frame attributes and data are cleared.
    /// @param frameType Frame type
    void setFrameType(unsigned int frameType);
    
    /// @brief Clears frame attributes and data
    void clear();
    
    /// @brief Finds an attribute and gives its value
    
    /// @param attrName Attribute name
    /// @param ok Pointer to a bool variable. If not @c CFNULL, the pointed 
    /// variable is set to @c true if the attribute was found; otherwise, it is
    /// set to @c false.
    /// @return Returns the attribute value (may be empty) or an empty string
    /// if the attribute was not found.
    std::string getAttribute(const std::string & attrName, bool * ok = CFNULL) const;
    
    /// @brief Checks if an attribute exists.
    
    /// @param attrName Attribute name
    /// @param emptyAllowed If @c false, empty attribute value is not allowed. 
    /// This means that if the attribute is found with an empty value, it is
    ///  considered as "not found".
    /// @return Returns @c true if the attribute was found, @c false if it was 
    /// not. 
    bool isAttributeSet(const std::string & attrName, bool emptyAllowed = false) const;
  };

//////////////////////////////////////////////////////////////////////////////

 }
}

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_BuilderParserFrameInfo_h