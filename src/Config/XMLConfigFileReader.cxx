// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <iostream>

#include "Common/StringOps.hh"
#include "Common/FilesystemException.hh"
#include "Config/ConfigArgs.hh"
#include "Config/XMLConfigFileReader.hh"
#include "Common/ParserException.hh"

#include "Config/ConverterTools.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

XMLConfigFileReader::XMLConfigFileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

XMLConfigFileReader::~XMLConfigFileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void XMLConfigFileReader::parseString(const std::string& xmlstr, ConfigArgs& args)
{
    // parse the file
    XMLResults pResults;
    XMLNode xMainNode = XMLNode::parseString(xmlstr.c_str(),"xml",&pResults);

    // display error message (if any)
    if (pResults.error != eXMLErrorNone)
    {
      std::string msg ("Error parsing xml string [") ; 
      msg += xmlstr;
	  msg += "]";
      throw Common::ParserException (FromHere(), msg);
    }

    process_xml_node(xMainNode,args);
}

//////////////////////////////////////////////////////////////////////////////

void XMLConfigFileReader::parse(const std::string& xcfcase, ConfigArgs& args)
{
    // parse the file
    XMLResults pResults;
    XMLNode xMainNode = XMLNode::parseFile(xcfcase.c_str(),"xml",&pResults);

    // display error message (if any)
    if (pResults.error != eXMLErrorNone)
    {
      std::string msg ("Error parsing file ") ; msg += xcfcase;
      throw Common::ParserException (FromHere(), msg);
    }

    args = ConverterTools::xCFcaseToConfigArgs(xMainNode);

}

//////////////////////////////////////////////////////////////////////////////

void XMLConfigFileReader::process_xml_node (XMLNode xNode, ConfigArgs& args, std::string nest)
{
 std::string node_name = xNode.getName();
 
 // add node name
 if (!nest.empty()) nest += node_name;
 // jump until the first XCFcase
 if (nest.empty() && node_name == "XCFcase") nest = node_name;
 // add separator
 if (!nest.empty()) nest += std::string(".");

 if (Common::StringOps::startsWith(nest,"XCFcase."))
 {
  std::string value;
  // each line is considered as a different text
  // so we need to set all of them in a unique string
  for(int i = 0 ; i < xNode.nText() ; i++)
  {
   value += xNode.getText(i);
    
   if(i < xNode.nText() - 1)
    value += '\n';
  }
   
  // if xNode does not have any child, it could be an option with an empty value
  if(!value.empty() || xNode.nChildNode() == 0)
   {
    std::string nested_label = nest.erase(nest.length() - 1, 1);
    // erase beginning "XCFcase."
    std::string nlabel = nested_label.erase(0,8);
    args[nested_label] = value;
   }
  }

  // call the method for each child
  for ( CFint nn = 0; nn < xNode.nChildNode(); ++nn)
  {
    process_xml_node ( xNode.getChildNode(nn), args, nest );
  }
  
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
