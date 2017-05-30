// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/DataHandleOutput.hh"
#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/TecWriter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<TecWriter,
               OutputFormatter,
               TecplotWriterModule,
               1>
tecWriterProvider("Tecplot");

//////////////////////////////////////////////////////////////////////////////

void TecWriter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("WriteSol","Command to wrtie the solution.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

TecWriter::TecWriter(const std::string& name)
  : OutputFormatter(name)
{
  addConfigOptionsTo(this);
  m_data.reset(new TecWriterData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);
  
  // default has been changed to be parallel
  m_writeSolutionStr = "ParWriteSolution";
  setParameter("WriteSol",&m_writeSolutionStr);
}
      
//////////////////////////////////////////////////////////////////////////////

TecWriter::~TecWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> TecWriter::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void TecWriter::setMethodImpl()
{
  OutputFormatter::setMethodImpl();

  setupCommandsAndStrategies();
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void TecWriter::unsetMethodImpl()
{
  cf_assert(m_unSetup.isNotNull());
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  OutputFormatter::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void TecWriter::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  OutputFormatter::configure(args);
  
  m_data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_data.getPtr(), args );
  
  // add configures to the TecWriterCom's

  configureCommand<TecWriterData,TecWriterComProvider>(args, m_setup, m_setupStr, m_data);

  configureCommand<TecWriterData,TecWriterComProvider>(args, m_unSetup, m_unSetupStr, m_data);

  configureCommand<TecWriterData,TecWriterComProvider>(args, m_writeSolution, m_writeSolutionStr, m_data);
  
  // automatically set m_appendRank=true if the write command name starts with "Par" 
  // (i.e. parallel output cases)
  if (StringOps::startsWith(m_writeSolutionStr, "Par")) {
    m_appendRank = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

void TecWriter::openImpl()
{
  CFAUTOTRACE;
  computeFullOutputName();
  m_data->setFilename(m_fullOutputName);
}

//////////////////////////////////////////////////////////////////////////////

void TecWriter::writeImpl()
{
  CFAUTOTRACE;
  cf_assert(isSetup());
  cf_assert(isConfigured());

  cf_assert(m_writeSolution.getPtr() != CFNULL);
  // set the nodal states at first
  getMethodData()->getCollaborator<SpaceMethod>()->extrapolateStatesToNodes();
  // write the solution file
  m_writeSolution->execute();
}

//////////////////////////////////////////////////////////////////////////////

void TecWriter::closeImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::string TecWriter::getFormatExtension() const
{
  CFAUTOTRACE;
  return std::string(".plt");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD
