// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Config/BuilderParserRules.hh"

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Config;

BuilderParserRules::BuilderParserRules(unsigned int errorType, unsigned int rootType, 
                                              const string & rootName)
{
 if(errorType == rootType)
  throw BuilderParserException(FromHere(), "Root and error types must be different.");
  
 if(rootName.empty())
  throw BuilderParserException(FromHere(), "Root name can not be empty.");
  
 this->m_errorType = errorType;
 this->m_rootType = rootType;
 this->m_frameNames[rootType] = rootName;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
BuilderParserRules::~BuilderParserRules() 
{
  // empty
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
bool BuilderParserRules::isValid(unsigned int type) const
{
 return this->m_frameNames.find(type)->first == type; 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
unsigned int BuilderParserRules::getRootType() const
{
 return this->m_rootType;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
unsigned int BuilderParserRules::getErrorType() const
{
 return this->m_errorType;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string BuilderParserRules::getRootName() const
{
 return this->m_frameNames.find(this->m_rootType)->second;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
string BuilderParserRules::getTypeName(unsigned int type) const
{
 if(this->isValid(type))
  return this->m_frameNames.find(type)->second;
 else
  return string();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
unsigned int BuilderParserRules::getFrameType(const string & name) const
{
 unsigned int foundType = this->m_errorType;
 map<unsigned int, string>::const_iterator it= this->m_frameNames.begin();
  
 while(it != this->m_frameNames.end() && foundType == this->m_errorType)
 {
  if(it->second == name)
   foundType = it->first;
  
  it++;
 }
  
 return foundType;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
BuilderParserFrameAttrs BuilderParserRules::getAttributes(unsigned int type, 
   bool * ok) const
{
 BuilderParserFrameAttrs attrs;
 bool valid = this->isValid(type);
 
 if(ok != CFNULL)
  *ok = valid;
  
 if(valid)
  attrs = this->m_frameAttributes.find(type)->second;
  
 return attrs;
}
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
bool BuilderParserRules::setTypeName(unsigned int type, const string & name)
{
 bool valid = !name.empty() && type != this->m_errorType;
   
 if(valid)
 {
  bool newItem = this->m_frameNames.find(type) == this->m_frameNames.end();
  this->m_frameNames[type] = name;
  
  if(newItem)
  {
   this->m_frameDataMand[type] = MAND_FORBIDDEN;
   this->m_frameAttributes[type] = BuilderParserFrameAttrs();
  }
 }
 return valid;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
bool BuilderParserRules::addAttribute(unsigned int type, const string & name, 
                                               BuilderParserMandType mandType)
{
 bool valid = this->check(type, name) && mandType != MAND_FORBIDDEN && 
   mandType != MAND_UNDEFINED;
  
 if(valid)
  this->m_frameAttributes[type][name] = mandType;
  
 return valid;
}
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
BuilderParserMandType BuilderParserRules::getAttrMandatoriness(unsigned int type, 
   const string & name) const
{
 BuilderParserMandType mandType = MAND_FORBIDDEN;
  
 if(this->check(type, name))
 {
  map<unsigned int, BuilderParserFrameAttrs>::const_iterator it;
  it = this->m_frameAttributes.find(type);
   
  if(it->first == type)
  {
   BuilderParserFrameAttrs attrs = it->second;
   BuilderParserFrameAttrs::const_iterator itAttrs = attrs.find(name);
   
   if(itAttrs->first == name)
    mandType = itAttrs->second;
  }
 }
  
 return mandType;
}
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
bool BuilderParserRules::setDataMandatoriness(unsigned int type, 
   BuilderParserMandType mandType)
{
 bool valid = this->isValid(type) && mandType != MAND_UNDEFINED;
  
 if(valid)
  this->m_frameDataMand[type] = mandType;
  
 return valid;
}
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
BuilderParserMandType BuilderParserRules::getDataMandatoriness(unsigned int type) const
{
 BuilderParserMandType mandType = MAND_UNDEFINED;
 map<unsigned int, BuilderParserMandType>::const_iterator it;
 
 it = this->m_frameDataMand.find(type);
 
 if(it != this->m_frameDataMand.end() && this->isValid(type))
  mandType = it->second;
 
 return mandType;
}

/****************************************************************************
 
                              PRIVATE METHODS
 
****************************************************************************/
 
bool BuilderParserRules::check(unsigned int type, const string & str) const
{
 return this->isValid(type) && !str.empty();
}
