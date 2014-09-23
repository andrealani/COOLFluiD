// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/FilterRHS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

FilterRHS::FilterRHS(const std::string& name) :
  NumericalStrategy(name)
{
  addConfigOptionsTo(this);

  m_maskIDs = std::vector<bool>();
  Config::ConfigObject::setParameter("maskIDs",&m_maskIDs);
}

//////////////////////////////////////////////////////////////////////////////

FilterRHS::~FilterRHS()
{
}

//////////////////////////////////////////////////////////////////////////////

void FilterRHS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<bool> >
    ("maskIDs","Flag telling if the current variable has to be filtered.");
}

//////////////////////////////////////////////////////////////////////////////

void FilterRHS::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
