// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MeshCreatorData.hh"
#include "Common/NoSuchValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void MeshCreatorData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("ScalingFactor","Mesh scaling factor");
  options.addConfigOption< std::string >("FileName","The file to open");
}

//////////////////////////////////////////////////////////////////////////////

MeshCreatorData::MeshCreatorData(Common::SafePtr<Method> owner)
 : MethodData(owner)
{
  addConfigOptionsTo(this);

  m_meshScalingFactor = 1.0;
  setParameter("ScalingFactor",&m_meshScalingFactor);

  m_fileName = "";
  setParameter("FileName",&m_fileName);
}

//////////////////////////////////////////////////////////////////////////////

MeshCreatorData::~MeshCreatorData()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreatorData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  MethodData::configure(args);

  if (m_fileName.empty())
  {
    throw Common::NoSuchValueException (FromHere(),"No filename specified!\n");
  }

  boost::filesystem::path fp (m_fileName);
  m_fileName = fp.string();
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreatorData::setup()
{
  MethodData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreatorData::unsetup()
{
  MethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

