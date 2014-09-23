// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <map>
#include <cstring> // for strcmp()

#include "Common/xmlParser.h"

#include "Config/BuilderParserFrameInfo.hh"
#include "Config/BuilderParserRules.hh"
#include "Config/BuilderParser.hh"

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Config;

string BuilderParser::m_errorString;

bool BuilderParser::buildFrame(const BuilderParserFrameInfo & frameInfo, 
                               const BuilderParserRules & rules,
                               string & frame)
{
 XMLNode rootNode;
 XMLNode typeNode;
 string typeName;
 BuilderParserMandType mand;
 bool hasChildren;
 BuilderParserFrameAttrs attributes;
 BuilderParserFrameAttrs::iterator it;
   
 // counter incremented each time we find a valid attribute. If, at the
 // end of the following loop, the counter is different from the number of 
 // provided attributes (i.e. the protocol defines that the type node needs 1 
 // attribute and 2 attributes are provided), this means that illegal 
 // (unknown) attributes were given, which is an error (protocol violation). 
 unsigned int attrCount = 0;
 unsigned int frameType = frameInfo.frameType;
 bool built = false;
 
 m_errorString.clear();
 
 // in order to avoid having many nested if-else block statements, each error 
 // is thrown as a exception and caught at the and of this function 
 try
 {
  if(!rules.isValid(frameType))
   throw "Unknown frame type";
  
  // create the root node
  rootNode = XMLNode::createXMLTopNode(rules.getRootName().c_str(), FALSE);
  
  // get the type name
  typeName = rules.getTypeName(frameType);
 
  // create the type node
  typeNode = rootNode.addChild(typeName.c_str()); 
   
  //
  // Check attributes
  //
 
  attributes = rules.getAttributes(frameType);
  it = attributes.begin();
   
  // if no attribute expected but at least one is provided
  if(attributes.empty() && !frameInfo.frameAttributes.empty())
   throw "Found attribute(s) although none was expected";
   
  while(it != attributes.end())
  {
   bool attrFound;
   string attrValue = frameInfo.getAttribute(it->first, &attrFound);
    
   // check if not found or the value is empty (only allowed if MAND_OPTIONAL)
   bool emptyValue = !attrFound || attrValue.empty();
    
   // if the value is empty but the attribute is mandatory
   if(emptyValue && it->second == MAND_MANDATORY)
    throw it->first + " attribute is mandatory and its value may not be empty";
 
   // (non-empty value) and (MAND_OPTIONAL or MAND_MANDATORY)
   if(!emptyValue) 
    typeNode.addAttribute(it->first.c_str(), attrValue.c_str());
    
   // no error, we increment the counter
   attrCount++;
    
   it++;
  }
   
  // if no error until now but the counter is different from the number of 
  // provided attributes
  if(attrCount != frameInfo.frameAttributes.size())
   throw "Found illegal attributes";
   
 
  //
  // Check data
  //
 
  mand = rules.getDataMandatoriness(frameType);
  hasChildren = frameInfo.frameData.nChildNode() > 0;
    
  if(mand == MAND_MANDATORY && !hasChildren)
   throw "Mandatory data not found";
    
  if(mand == MAND_FORBIDDEN && hasChildren)
   throw "Data forbidden for this type";
    
  if((mand == MAND_MANDATORY || mand == MAND_OPTIONAL) && hasChildren)
   typeNode.addChild(frameInfo.frameData);
 
  built = true;
  frame = rootNode.createXMLString(); 
 }
 catch(const char * e)
 {
  built = false;
  m_errorString = e;
 }
 
 return built;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool BuilderParser::parseFrame(const string & xmlString,
                               const BuilderParserRules & rules,
                               BuilderParserFrameInfo & frameInfo)
{
 XMLNode node;
 XMLResults error;
 XMLNode typeNode;
 unsigned int frameType;
 BuilderParserFrameAttrs attrs;
 BuilderParserFrameAttrs::iterator it;
 BuilderParserMandType mand;
 bool hasChildren;
  
 bool valid = false;
 node = XMLNode::parseString(xmlString.c_str(), rules.getRootName().c_str(), &error);
 
 // attrsCount has the same utility as in BuilderParser::buildFrame()
 int attrCount = 0;

 m_errorString.clear();
 
 // in order to avoid having many nested if-else block statements, each error 
 // is thrown as a exception and caught at the and of this function 
 try
 {
  //
  // Check root and type nodes
  //
  
  if(error.error != eXMLErrorNone)
   throw XMLNode::getError(error.error);
  
  if(strcmp(node.getName(), rules.getRootName().c_str()) != 0)
   throw "Root node was not found";
   
  if(node.nChildNode() != 1)
   throw "Root node must have exactly one child node";
   
  typeNode = node.getChildNode(0);
  frameType = rules.getFrameType(typeNode.getName());
    
  if(frameType == rules.getErrorType())
   throw "Unknown frame type";
 
  frameInfo.setFrameType(frameType); 
  attrs = rules.getAttributes(frameType);
  it = attrs.begin();
  
  //
  // Check attributes
  //
  
  while(it != attrs.end())
  {
   const char * value = typeNode.getAttribute(it->first.c_str());
   
   // if attribute is mandatory but has an empty value (or was not found)
   if(it->second == MAND_MANDATORY && (value == NULL || strlen(value) == 0))
    throw string(it->first + " mandatory attribute is missing or has an empty value.").c_str();
   
   if(value != NULL && strlen(value) != 0)
   {
    frameInfo.frameAttributes[it->first] = value;
    attrCount++;
   }
   
   it++;
  }
  
  if(attrCount != typeNode.nAttribute() && m_errorString.empty())
   throw "Found illegal attributes.";
  
  //
  // Check data
  //
  
  mand = rules.getDataMandatoriness(frameType);
  hasChildren = typeNode.nChildNode() > 0;
    
  if(mand == MAND_MANDATORY && !hasChildren)
   throw "Mandatory data not found";
      
  if(mand == MAND_FORBIDDEN && hasChildren)
   throw "Data forbidden for this type";
      
  if((mand == MAND_MANDATORY || mand == MAND_OPTIONAL) && hasChildren)
  {
   // clear frame info data and create a new node
   frameInfo.frameData.deleteNodeContent();
   frameInfo.frameData = XMLNode::createXMLTopNode("", FALSE);

   // append all children to the frame data node
   for(int i = 0 ; i < typeNode.nChildNode() ; i++)
    frameInfo.frameData.addChild(typeNode.getChildNode(i).deepCopy());
  }
      
  valid = true;
 }
 catch(const char * e)
 {
  valid = false;
  frameInfo.setFrameType(rules.getErrorType());
  m_errorString = e;
 }
 
 return valid;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void BuilderParser::createXMLNode(XMLNode & xmlNode)
{
 xmlNode.deleteNodeContent(); 
 xmlNode = XMLNode::createXMLTopNode("xml", TRUE); 
  
 xmlNode.addAttribute("version","1.0");
 xmlNode.addAttribute("encoding","UTF-8");
 xmlNode.addAttribute("standalone","yes");
}

/****************************************************************************

                              PRIVATE FUNCTIONS

****************************************************************************/

string BuilderParser::getErrorString()
{
 return m_errorString;
}