// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullOutputFormatter.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullOutputFormatter,
               OutputFormatter,
               FrameworkLib,
               1>
nullOutputFormatterProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullOutputFormatter::NullOutputFormatter(const std::string& name)
  : OutputFormatter(name)
{
  m_filenameStr = "null";
}

//////////////////////////////////////////////////////////////////////////////

NullOutputFormatter::~NullOutputFormatter()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullOutputFormatter::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullOutputFormatter::setMethodImpl()
{
  OutputFormatter::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullOutputFormatter::unsetMethodImpl()
{
  OutputFormatter::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullOutputFormatter::openImpl()
{
  CFLog(VERBOSE,"NullOutputFormatter::open() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullOutputFormatter::writeImpl()
{
  CFLog(VERBOSE,"NullOutputFormatter::write() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullOutputFormatter::closeImpl()
{
  CFLog(VERBOSE,"NullOutputFormatter::close() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

std::string NullOutputFormatter::getFormatExtension() const
{
  CFLog(VERBOSE,"NullOutputFormatter::getFormatExtension() called!" << "\n");
  return std::string(".null");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
