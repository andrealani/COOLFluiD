// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CommandGroup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void CommandGroup::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("groupedComs","Name of the NumericalCommand belonging to this CommandGroup.");
   options.addConfigOption< std::vector<std::string> >("groupedTRS","Names of the TRS's belonging to this CommandGroup.");
}

//////////////////////////////////////////////////////////////////////////////

CommandGroup::CommandGroup(const std::string& name)
   : ConfigObject(name)
{
   addConfigOptionsTo(this);
  _trsNames = std::vector<std::string>();
   setParameter("groupedTRS",&_trsNames);


  _comsNames = std::vector<std::string>();
   setParameter("groupedComs",&_comsNames);


}

//////////////////////////////////////////////////////////////////////////////

CommandGroup::~CommandGroup()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
