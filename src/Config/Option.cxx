// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>

#include "Config/Option.hh"
#include "Config/OptionValidation.hh"
#include "Config/OptionValidationException.hh"
#include "Config/BadMatchException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

const string Option::DELIMITER = "--";

//////////////////////////////////////////////////////////////////////////////

Option::Option(const std::string& label,
               const std::string& description) :
  m_label(label),
  m_description(description),
  m_owner_hname(),
  m_dynamic(false),
  m_basic(false)
{
}

//////////////////////////////////////////////////////////////////////////////

Option::Option(const Option& other) :
  m_label(other.m_label),
  m_description(other.m_description),
  m_owner_hname(),                     // might have different parent
  m_dynamic(other.m_dynamic),
  m_basic(other.m_basic)
{
}

//////////////////////////////////////////////////////////////////////////////

Option::~Option()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string Option::getLabel() const
{
  return m_label;
}

//////////////////////////////////////////////////////////////////////////////

void Option::setLabel(const std::string& label)
{
  m_label = label;
}

//////////////////////////////////////////////////////////////////////////////

std::string Option::getNestedLabel() const
{
  return m_owner_hname + m_label;
}

//////////////////////////////////////////////////////////////////////////////

void Option::setHierarchyOwner(const std::string& owner)
{
  m_owner_hname = owner;
}

//////////////////////////////////////////////////////////////////////////////

std::string Option::getHierarchyOwner() const
{
  return m_owner_hname;
}

//////////////////////////////////////////////////////////////////////////////

std::string Option::description() const
{
  return m_description;
}

//////////////////////////////////////////////////////////////////////////////

std::string Option::helpDescription() const
{
  std::ostringstream oss;
  oss << getLabel() << "\t";
  oss << description();
  oss << "\n";
  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

void Option::onProcessed()
{
}

//////////////////////////////////////////////////////////////////////////////

void Option::onUpdate()
{
}

//////////////////////////////////////////////////////////////////////////////

void Option::onComplete(const OptionList& opts)
{
}

//////////////////////////////////////////////////////////////////////////////

bool Option::matchNestedLabel(const std::string& label) const
{
  return label == getNestedLabel();
}

//////////////////////////////////////////////////////////////////////////////

void Option::addLinkedOption(Option* const opt)
{
  m_linked_opts.push_back(opt);
}

//////////////////////////////////////////////////////////////////////////////

vector<Option*>& Option::getLinkedOptions()
{
  return m_linked_opts;
}

//////////////////////////////////////////////////////////////////////////////

bool Option::matchesDelimiter(const std::string& arg)
{
  return arg.substr(0, Option::DELIMITER.length()) == Option::DELIMITER;
}

//////////////////////////////////////////////////////////////////////////////

void Option::checkArgument(const string& arg)
{
  if (matchesDelimiter(arg))
  {
    string errmsg = "argument '" + arg
                  + "' for option " + m_label
                  + " appears to be a tag";
    throw BadMatchException (FromHere(),errmsg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Option::setLinkedOptions(const std::string& value)
{
  for (std::vector<Option*>::iterator it = getLinkedOptions().begin(), stop = getLinkedOptions().end();
       it != stop; ++it)
  {
    Option* opt = *it;
    opt->setValue(value);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Option::clearValidations()
{
  for (std::vector<OptionValidation*>::iterator
         it   = m_validations.begin(),
         stop = m_validations.end();
         it != stop; ++it)
  {
    deletePtr(*it);
  }
  m_validations.clear();
}

//////////////////////////////////////////////////////////////////////////////

void Option::validate ()
{
  for (std::vector<OptionValidation*>::iterator
         it   = m_validations.begin(),
         stop = m_validations.end();
         it != stop; ++it)
  {
    OptionValidation* val = *it;
    if ( !val->isValid ())
    {
      std::ostringstream oss;
      oss << "Option [" << getLabel()
          << " with value [" << getValueAsString()
          <<  "] is not valid because " << val->getReason();
      throw OptionValidationException (FromHere(),oss.str());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Option::addValidation(OptionValidation* val)
{
  m_validations.push_back(val);
}

//////////////////////////////////////////////////////////////////////////////

void Option::makeOptionDynamic()
{
  m_dynamic = true;
}

//////////////////////////////////////////////////////////////////////////////

void Option::makeOptionBasic()
{
  m_basic = true;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD
