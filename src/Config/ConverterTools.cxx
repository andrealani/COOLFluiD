// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>

#include <string>
#include <vector>

#include "Common/StringOps.hh"
#include "Common/xmlParser.h"

#include "Config/ConfigArgs.hh"
#include "Common/ParserException.hh"

#include "Config/ConverterTools.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace std;

std::string ConverterTools::configArgsToXml(const ConfigArgs & args)
{
 return ConverterTools::configArgsToXmlNode(args).createXMLString(); 
}

//////////////////////////////////////////////////////////////////////////////

std::string ConverterTools::configArgsToCFcase(const ConfigArgs & args)
{
 std::string cfcase;
 ConfigArgs::const_iterator itArgs = args.begin();
 
 for( ; itArgs != args.end() ; itArgs++)
  cfcase += itArgs->first + " = " + itArgs->second + '\n';
 
 return cfcase;
}

//////////////////////////////////////////////////////////////////////////////

XMLNode ConverterTools::configArgsToXCFcase(const ConfigArgs & args)
{
 XMLNode rootNode = XMLNode::createXMLTopNode("xml", TRUE);
 XMLNode xcfcaseNode;
 
 rootNode.addAttribute("version","1.0");
 rootNode.addAttribute("encoding","UTF-8");
 rootNode.addAttribute("standalone","yes");
 
 xcfcaseNode = rootNode.addChild("XCFcase");
 
 xcfcaseNode.addChild(ConverterTools::configArgsToXmlNode(args));
 
 return rootNode;
}

//////////////////////////////////////////////////////////////////////////////

ConfigArgs ConverterTools::xmlToConfigArgs(const std::string & xmlStr)
{
 ConfigArgs args;
 XMLNode xMainNode = ConverterTools::xmlToXMLNode(xmlStr);
 
 ConverterTools::xmlNodeToArgs(xMainNode, args);
 
 return args;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConverterTools::xmlToCFcase(const std::string & xmlStr)
{
 return ConverterTools::configArgsToCFcase(ConverterTools::xmlToConfigArgs(xmlStr));
}

//////////////////////////////////////////////////////////////////////////////

XMLNode ConverterTools::xmlToXCFcase(const std::string & xmlStr)
{
 XMLNode xMainNode = ConverterTools::xmlToXMLNode(xmlStr);
 XMLNode rootNode = XMLNode::createXMLTopNode("xml", TRUE);
 XMLNode xcfcaseNode;
 
 rootNode.addAttribute("version","1.0");
 rootNode.addAttribute("encoding","UTF-8");
 rootNode.addAttribute("standalone","yes");
 
 xcfcaseNode = rootNode.addChild("XCFcase");
 
 xcfcaseNode.addChild(xMainNode);
 
 return rootNode; 
}

//////////////////////////////////////////////////////////////////////////////

XMLNode ConverterTools::createNode(const XMLNode & root, const vector<string> & path)
{
 vector<string>::const_iterator itPath = path.begin();
 XMLNode currentNode = root;
 
 for( ; itPath != path.end() ; itPath++)
 {
  string nodeName = *itPath;
   
  XMLNode node = currentNode.getChildNode(nodeName.c_str(), 0); 
  
  // if the node does not exist we create it
  if (node.isEmpty()) 
   node = currentNode.addChild(nodeName.c_str());
  
  currentNode = node;
 }
 return currentNode; 
}

//////////////////////////////////////////////////////////////////////////////

XMLNode ConverterTools::configArgsToXmlNode(const ConfigArgs & args)
{
 ConfigArgs::const_iterator itArgs = args.begin();
 XMLNode rootNode = XMLNode::createXMLTopNode("", FALSE);

 for( ; itArgs != args.end() ; itArgs++)
 {
  vector<string> path = StringOps::getWords(itArgs->first, '.');
  XMLNode node = ConverterTools::createNode(rootNode, path);
  
  node.addText(itArgs->second.c_str());
 }
 
 return rootNode; 
}

//////////////////////////////////////////////////////////////////////////////

XMLNode ConverterTools::xmlToXMLNode(const std::string & xmlStr)
{
 XMLResults pResults;
 XMLNode xMainNode = XMLNode::parseString(xmlStr.c_str(), "", &pResults);
 
 // display error message (if any)
 if (pResults.error != eXMLErrorNone)
 {
  std::string msg ("Error parsing xml string [") ; 
  msg += xmlStr;
  msg += "]: ";
  msg += XMLNode::getError(pResults.error);
  throw Common::ParserException (FromHere(), msg);
 }
 
 return xMainNode; 
}

//////////////////////////////////////////////////////////////////////////////

void ConverterTools::xmlNodeToArgs(const XMLNode & xNode, ConfigArgs& args, std::string nest)
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
  xmlNodeToArgs( xNode.getChildNode(nn), args, nest );
 }
  
}

//////////////////////////////////////////////////////////////////////////////

std::string ConverterTools::xCFcaseToXml(const XMLNode & xNode)
{
 XMLNode xcfcaseNode = xNode.getChildNode("XCFcase", 0);
 char * xmlCharPtr;
 std::string returnString;
 
 if(!xcfcaseNode.isEmpty() && xcfcaseNode.nChildNode() != 0)
  xmlCharPtr = xcfcaseNode.getChildNode(0).createXMLString();
 
 xmlCharPtr = xNode.createXMLString();
 
 returnString = xmlCharPtr;
 delete xmlCharPtr;
 return returnString;
}

//////////////////////////////////////////////////////////////////////////////

ConfigArgs ConverterTools::xCFcaseToConfigArgs(const XMLNode & xNode)
{
 ConfigArgs args;
 
 if(!xNode.isEmpty())
  ConverterTools::xmlNodeToArgs(xNode, args);
 
 return args;
}
