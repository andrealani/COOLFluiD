// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_BuilderParser_h
#define COOLFluiD_Config_BuilderParser_h

/////////////////////////////////////////////////////////////////////////////

#include "Config/Config.hh"
#include <string>

class XMLNode;

namespace COOLFluiD {
 
 namespace Config {

/////////////////////////////////////////////////////////////////////////////

  class BuilderParserFrameInfo;
  
  class BuilderParserRules;
  
  /// @brief Builds and parses a frame according to a given protocol.
  
  /// @author Quentin Gasper.
  
  class Config_API BuilderParser
  {
   public:
    
    /// @brief Builds a frame
    
    /// @param frameInfo Frame contents.
    /// @param rules Protocol the built frame has to respect.
    /// @param frame String were the built frame will be stored. The string is 
    /// cleared at the beginning of the function execution.
    /// @return Returns @c true if the frame was successfully built. Otherwise
    /// returns @c false. In the latter case, error message can be obtained
    /// with @c #getErrorString() and @c frame is empty.
    static bool buildFrame(const BuilderParserFrameInfo & frameInfo,
                           const BuilderParserRules & rules,
                           std::string & frame);
    
    /// @brief Parses a frame
    
    /// @param xmlString String to parse.
    /// @param rules Protocol the built frame has to respect
    /// @param frameInfo Structure where frame contents will be stored. The 
    /// object is cleared at the beginning of the function execution.
    /// @return Returns @c true if the frame was successfully parsed. Otherwise
    /// returns @c false. In the latter case, error message can be obtained
    /// with @c #getErrorString() and @c frameInfo is empty.
    static bool parseFrame(const std::string & xmlString,
                           const BuilderParserRules & rules,
                           BuilderParserFrameInfo & frameInfo);
    
    /// @brief Gives the last error message.
    
    /// @return Returns the last error message.
    static std::string getErrorString();
    
    /// @brief Creates a XML node.
    
    /// The built node is @code &lt;?xml version="1.0" encoding="UTF-8" 
    /// standalone="yes" ?&gt; @endcode 
    /// @param xmlNode Object were the new node will be stored. Object 
    /// contents are deleted.
    static void createXMLNode(XMLNode & xmlNode);
    
   private:
    
    /// @brief Last error message.
    static std::string m_errorString;
    
  }; // class BuilderParser

/////////////////////////////////////////////////////////////////////////////

 } // namespace Config
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_BuilderParser_h