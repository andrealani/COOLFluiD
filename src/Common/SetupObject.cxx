// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/SetupObject.hh"
#include "Common/SetupException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

SetupObject::SetupObject() : m_setup(false)
{
}

//////////////////////////////////////////////////////////////////////////////

SetupObject::~SetupObject()
{
}

//////////////////////////////////////////////////////////////////////////////

bool SetupObject::isSetup()
{
  return m_setup;
}

//////////////////////////////////////////////////////////////////////////////

void SetupObject::setup()
{
  if (m_setup)
    throw SetupException (FromHere(),"Trying to setup object again without calling unsetup first.");
  m_setup = true;
}

//////////////////////////////////////////////////////////////////////////////

void SetupObject::unsetup()
{
  if (!m_setup)
    throw SetupException (FromHere(),"Trying to unsetup object without calling setup first.");
  m_setup = false;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

