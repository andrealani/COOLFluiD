// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/OutputFormatterData.hh"
#include "Framework/DataHandleOutput.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void OutputFormatterData::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

OutputFormatterData::OutputFormatterData(Common::SafePtr<Method> owner)
 : MethodData(owner)
{
  addConfigOptionsTo(this);

  m_datah_out = new DataHandleOutput("DataHandleOutput");
}

//////////////////////////////////////////////////////////////////////////////

OutputFormatterData::~OutputFormatterData()
{
  deletePtr(m_datah_out);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<DataHandleOutput> OutputFormatterData::getDataHOutput()
{
  return m_datah_out;
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatterData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  MethodData::configure(args);

  m_datah_out->setParentNamespace(getNamespace());
  m_strategies.push_back ( m_datah_out );
  configureNested ( m_datah_out, args );
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatterData::setup()
{
  MethodData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatterData::unsetup()
{
  MethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

