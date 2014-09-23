// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_ConverterTools_h
#define COOLFluiD_Config_ConverterTools_h

//////////////////////////////////////////////////////////////////////////////

#include "Config/Config.hh"
#include "Common/xmlParser.h"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
 namespace Config {

//////////////////////////////////////////////////////////////////////////////

  class ConfigArgs;
  
  /// @brief This class provides tools for conversions between ConfigArgs, 
  /// CFcase, XML and XCFcase formats.
  
  /// @author Quentin Gasper
  
  class Config_API ConverterTools
  {   
   public:
    /// @brief Converts a ConfigArgs object to an XML string
    
    /// @param args ConfigArgs to convert
    /// @return Returns the built XML string
    static std::string configArgsToXml(const ConfigArgs & args);
    
    /// @brief Converts a ConfigArgs object to an CFcase string
    
    /// @param args ConfigArgs to convert
    /// @return Returns the built CFcase string
    static std::string configArgsToCFcase(const ConfigArgs & args);
    
    /// @brief Converts a ConfigArgs object to an XML string
    
    /// @param args ConfigArgs to convert
    /// @return Returns the built XML string
    static XMLNode configArgsToXCFcase(const ConfigArgs & args);
    
    /// @brief Converts an XML string to a ConfigArgs object
    
    /// @param xmlStr XML string to convert
    /// @return Returns the built ConfigArgs object
    static ConfigArgs xmlToConfigArgs(const std::string & xmlStr);
    
    /// @brief Converts an XML string to a CFcase string
    
    /// @param xmlStr XML string to convert
    /// @return Returns the built CFcase string
    static std::string xmlToCFcase(const std::string & xmlStr);
    
    /// @brief Converts an XML string to an XML node
    
    /// @param xmlStr XML string to convert
    /// @return Returns the built XML node
    static XMLNode xmlToXCFcase(const std::string & xmlStr);
    
    /// @brief Converts an XMLNode object to an XML string
    
    /// If the provided node has a child named "XCFcase", this method returns 
    /// the string representation of the first child of this "XCFcase" node. 
    /// If the provided node does not have a such child, this method returns
    /// its string representation.
    /// If the node is empty, this method returns an empty string.
    /// @param xNode XML node to convert
    /// @return Returns the built XML string.
    static std::string xCFcaseToXml(const XMLNode & xNode);
    
    /// @brief Converts an XMLNode object to a ConfigArgs
    
    /// @param xNode XML node to convert
    /// @return Returns the built ConfigArgs object
    static ConfigArgs xCFcaseToConfigArgs(const XMLNode & xNode);
    
    /// @brief Converts an XMLNode object to a CFcase string

    /// @param xNode XML node to convert
    /// @return Returns the built CFcase string
    static ConfigArgs cfcaseToConfigArgs(const std::string & cfcase);
   
   private:
    /// @brief Creates a node.
    
    /// The path indicate the new node path; its last element is the new node
    /// name. The path is relative to the root, which means that the first
    /// element of the vector is a child of the root. If a node in the path
    /// does not exist, it is created.
    /// @param root Root node to which the path is relative
    /// @param path Path of the new node
    /// @return Returns the created node, or the root node if the path is empty.
    static XMLNode createNode(const XMLNode & root, const std::vector<std::string> & path);
    
    /// @brief Converts a ConfigArgs object to XML.
    
    /// @param args ConfigArgs object to convert
    /// @return Returns the root node of the built XML tree.
    static XMLNode configArgsToXmlNode(const ConfigArgs & args);
    
    /// @brief Converts an XML string to an XMLNode object.
    
    /// @param xmlStr XML string to convert.
    /// @return Returns the root node of the built XML tree.
    /// @throws ParserException If any error occurs during the string parsing.
    static XMLNode xmlToXMLNode(const std::string & xmlStr);
    
    /// @brief Converts an XML node to a ConfigArgs object
    
    /// The code of this function comes from XMLConfigFileReader::process_xml_node()
    /// @param xNode Node to convert.
    /// @param args ConfigArgs where the data will be stored.
    /// @param nest Contains parent node names, separated by '.'.
    static void xmlNodeToArgs(const XMLNode & xNode, ConfigArgs& args, std::string nest = "");
 
  }; // class ConverterTools
 
////////////////////////////////////////////////////////////////////////////

 } // namespace Config
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_ConverterTools_h