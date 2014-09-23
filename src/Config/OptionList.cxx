// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>

#include "Config/OptionList.hh"
#include "Config/BadMatchException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

OptionList::OptionList() : m_strictArgs(false) {}

//////////////////////////////////////////////////////////////////////////////

OptionList::OptionList(const OptionList& obj)
{
  operator=(obj);
}

//////////////////////////////////////////////////////////////////////////////

const OptionList& OptionList::operator=(const OptionList& obj)
{
  deleteOptions();
  std::vector<Option*>::const_iterator itr = obj.m_vec_options.begin();
  for(; itr != obj.m_vec_options.end(); ++itr) {
    cf_assert(*itr != CFNULL);
    m_vec_options.push_back((*itr)->clone());
  }
  m_strictArgs = obj.m_strictArgs;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::deleteOptions()
{
  vector<Option*>::iterator it = m_vec_options.begin();
  vector<Option*>::iterator stop = m_vec_options.end();
  for (; it != stop; ++it) {
    Option* opt = *it;
    deletePtr(opt);
  }
  m_vec_options.clear();
}

//////////////////////////////////////////////////////////////////////////////

OptionList::~OptionList()
{
  deleteOptions();
}

//////////////////////////////////////////////////////////////////////////////

std::pair<ConfigArgs,ConfigArgs>
OptionList::processArgs(ConfigArgs& args, bool dynamic)
{

  std::pair< ConfigArgs, ConfigArgs > result;

  std::vector<ConfigKey> used_args;

  // loop over all args supplied
  ConfigArgs::const_iterator begin = args.begin();
  ConfigArgs::const_iterator end   = args.end();
  for (ConfigArgs::const_iterator argIt = begin; argIt != end; ++argIt)
  {
    bool consumed = false;

    std::string label = argIt->first;
    std::string value = argIt->second;

    // loop over all options for matches
    vector<Option*>::iterator start = m_vec_options.begin();
    vector<Option*>::iterator stop  = m_vec_options.end();

    for (vector<Option*>::iterator optIt = start; optIt != stop; ++optIt)
    {
      Option* opt = *optIt;

      // skip if option is not dynamic and the arguments are
      if (dynamic && !opt->isDynamic()) continue;

      if (opt->matchNestedLabel(label))
      {
        vector<std::string> matches = findMatches(label);

        if (matches.size() > 1) { complainAboutMatches(matches, label); }

        consumed = true;

        used_args.push_back(label);

        opt->setValue(value);

        if (dynamic && opt->isDynamic())
          opt->onUpdate();
        else
          opt->onProcessed();

      } // option mathced

    } // loop all options

    // the argument was consumed ?
    if (consumed)
      result.first[label] = value;
    else
      result.second[label] = value;

  } // loop arguments

  args.consume(used_args);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::processCommandLine (int argc, char** argv)
{
  // skip the first argument, since that's the name of the program:
  int i = 1;

  while (i < argc)
  {
    if (Option::matchesDelimiter(argv[i]))
    {
      // option without delimiter
      std::string arg = std::string(std::string(argv[i]), Option::DELIMITER.length());

      vector<Option*>::iterator start = m_vec_options.begin();
      vector<Option*>::iterator stop = m_vec_options.end();
      bool consumed = false;

      for (vector<Option*>::iterator it = start; it != stop; ++it)
      {
        Option* opt = *it;

        if (opt->matchNestedLabel(arg))
        {
          vector<std::string> matches = findMatches(arg);

          if (matches.size() > 1) { complainAboutMatches(matches,arg); }

          consumed = true;

          i += opt->setCommandLineValues(arg, i, argc, argv);

          // not doing anything with the return value here.
          opt->onProcessed();
        }
      }

      if (!consumed)
      {
        string argstr = string(argv[i]);
        // maybe it was a cry for help
        if (argstr == "--help" ||
            argstr == "--h"    ||
            argstr == "-h"     ||
            argstr == "--?"    ||
            argstr == "-?"     ||
            argstr == "?")
        {
            std::cout << writeUsage();
            return;
        }
        if (argstr == "--help-config")
        {
            std::cout << writeConfigUsage();
            return;
        }
        // the argument was not consumed -- did it look like an option?
        else if (getStrictArgs())
        {
            string msg("no option '" + string(argv[i]) + string("'"));
            throw BadMatchException (FromHere(),msg);
        }
        else { /* argument not consumed */ }

//        m_unprocessed_args.push_back(argv[i]);
        i++;
      }
    }
    else {
      // arg was not consumed, and does not look like an option
//      m_unprocessed_args.push_back(argv[i]);
      i++;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::string OptionList::writeUsage() const
{
  ostringstream oss;

  // length of the longest name, for alignment
  size_t maxlen = 0;

  vector<Option*>::const_iterator it = m_vec_options.begin();
  vector<Option*>::const_iterator stop = m_vec_options.end();

  for (; it != stop; ++it)
  {
    Option* opt = *it;
    maxlen = max(maxlen, opt->getLabel().length());
  }

  // give it some more space (literally)
  maxlen += 4;

  // reset to head.
  it = m_vec_options.begin();
  for (; it != stop; ++it)
  {
    Option* opt = *it;
    std::string nm = opt->getLabel();
    oss << "    ";
    oss << nm;
    size_t len = nm.length();
    for (size_t i = len; i < maxlen; ++i) {  oss << " "; }
    oss << opt->helpDescription();
    oss << "\n";
  }
  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

std::string OptionList::writeConfigUsage() const
{
  ostringstream oss;

  vector<Option*>::const_iterator it = m_vec_options.begin();
  vector<Option*>::const_iterator stop = m_vec_options.end();

  for (; it != stop; ++it)
  {
    Option* opt = *it;
    oss << "# ";
    // should probably insert # after newlines.
    oss << opt->description();
    oss << "\n";

    oss << opt->getLabel() << " = " << opt->getDefaultValueAsString() << "\n";
    oss << "\n";
  }
  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

Option* OptionList::addOptionPtr(Option* opt)
{
  checkNameDoesNotExist(opt->getLabel());
  m_vec_options.push_back(opt);
  return opt;
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::linkOptions(const std::string& from, const std::string& to)
{
  Option* fromopt = findFirstByName(from);
  Option* toopt   = findFirstByName(to);
  if ((fromopt != CFNULL) && (toopt != CFNULL)) { toopt->addLinkedOption(fromopt); }
  else
  {
    std::string msg = "Tryingo to link inexisting options: [" + from + "] or [" + to + "]";
    throw ConfigOptionException (FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::setHierarchyOptions(const std::string& nest)
{
  vector<Option*>::iterator it = m_vec_options.begin();
  vector<Option*>::iterator stop =  m_vec_options.end();
  for (; it != stop; ++it) {
    Option* opt = *it;
    opt->setHierarchyOwner(nest);
  }
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::prependToName(const std::string& pre)
{
  vector<Option*>::iterator it = m_vec_options.begin();
  vector<Option*>::iterator stop =  m_vec_options.end();
  for (; it != stop; ++it) {
    Option* opt = *it;
    if(!Common::StringOps::startsWith(opt->getLabel(),pre)) {
      opt->setLabel(pre + opt->getLabel());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<Option*> OptionList::findByName(const std::string& name) const
{
  vector<Option*> matches;

  vector<Option*>::const_iterator it = m_vec_options.begin();
  vector<Option*>::const_iterator stop = m_vec_options.end();

  for (; it != stop; ++it)
  {
    Option* opt = *it;
    if (opt->getLabel() == name)
    {
      matches.push_back(opt);
    }
  }
  return matches;
}

//////////////////////////////////////////////////////////////////////////////

Option* OptionList::findFirstByName(const std::string& name) const
{
  vector<Option*> matches = findByName(name);
  return matches.size() > 0 ? matches[0] : CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::checkNameDoesNotExist(const std::string& name) const
{
  vector<Option*> matches = findByName(name);
  if (matches.size() > 0)
  {
     std::string errmsg ("Error: a config option with the name : ");
     errmsg += name;
     errmsg += " already matches ";
     errmsg += matches.size();
     errmsg += " other option(s) ";
    throw DuplicateNameException (FromHere(),errmsg);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> OptionList::findMatches(const std::string& label) const
{
  vector<std::string> matches;

  vector<Option*>::const_iterator it = m_vec_options.begin();
  vector<Option*>::const_iterator stop = m_vec_options.end();
  for (; it != stop; ++it) {
    Option* opt = *it;
    if (opt->matchNestedLabel(label)) {
      matches.push_back(opt->getLabel());
    }
  }
  return matches;
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::complainAboutMatches(const vector<std::string>& matches, const std::string& tag) const
{
  if (matches.size() > 1) {
    std::string errmsg = tag + " is ambiguous. ";
    errmsg += "It could match ";
    if (matches.size() == 2) {
      errmsg += "either " + matches[0] + " or " + matches[1];
    }
    else {
      CFuint nMatches = matches.size();
      for (CFuint mi = 0; mi < nMatches - 2; ++mi) {
        errmsg += matches[mi] + ", ";
      }
      errmsg += matches[nMatches - 2] + " or " + matches[nMatches - 1];
    }
    throw BadMatchException (FromHere(),errmsg);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::string OptionList::debugInfo() const
{
  ostringstream oss;

  vector<Option*>::const_iterator it = m_vec_options.begin();
  vector<Option*>::const_iterator stop = m_vec_options.end();

  oss << "nb options : " << m_vec_options.size() << "\n";
  for (; it != stop; ++it)
  {
    Option* opt = *it;
    oss << opt->debugInfo();
    oss << "\n";
  }
  return oss.str();
}


//////////////////////////////////////////////////////////////////////////////

void OptionList::doComplete()
{
  vector<Option*>::iterator oit = m_vec_options.begin();
  vector<Option*>::iterator ostop = m_vec_options.end();
  for (; oit != ostop; ++oit)
  {
      Option* opt = *oit;
      opt->onComplete(*this);
  }
}

//////////////////////////////////////////////////////////////////////////////

void OptionList::doValidation()
{
  vector<Option*>::iterator oit = m_vec_options.begin();
  vector<Option*>::iterator ostop = m_vec_options.end();
  for (; oit != ostop; ++oit)
  {
      Option* opt = *oit;
      opt->validate();
  }
}

//////////////////////////////////////////////////////////////////////////////

std::string OptionList::getOptionsXML() const
{
  std::string result;
  vector<Option*>::const_iterator it   = m_vec_options.begin();
  vector<Option*>::const_iterator stop = m_vec_options.end();

  for (; it != stop; ++it)
    result += (*it)->getXML();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
