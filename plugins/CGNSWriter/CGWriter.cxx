// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/DataHandleOutput.hh"
#include "CGNSWriter/CGNSWriter.hh"
#include "CGNSWriter/CGWriter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CGWriter,
               OutputFormatter,
               CGNSWriterModule,
               1>
CGWriterProvider("CGNS");

//////////////////////////////////////////////////////////////////////////////

void CGWriter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("WriteSol","Command to wrtie the solution.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

CGWriter::CGWriter(const std::string& name)
  : OutputFormatter(name)
{
  addConfigOptionsTo(this);
  m_data.reset(new CGWriterData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);
  
  // default has been changed to be parallel
  m_writeSolutionStr = "ParWriteSolution";
  setParameter("WriteSol",&m_writeSolutionStr);
}
      
//////////////////////////////////////////////////////////////////////////////

CGWriter::~CGWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> CGWriter::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void CGWriter::setMethodImpl()
{
  OutputFormatter::setMethodImpl();

  setupCommandsAndStrategies();
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void CGWriter::unsetMethodImpl()
{
  cf_assert(m_unSetup.isNotNull());
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  OutputFormatter::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void CGWriter::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  OutputFormatter::configure(args);
  
  m_data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_data.getPtr(), args );
  
  // add configures to the CGWriterCom's

  configureCommand<CGWriterData,CGWriterComProvider>(args, m_setup, m_setupStr, m_data);

  configureCommand<CGWriterData,CGWriterComProvider>(args, m_unSetup, m_unSetupStr, m_data);

  configureCommand<CGWriterData,CGWriterComProvider>(args, m_writeSolution, m_writeSolutionStr, m_data);
  
  // automatically set m_appendRank=true if the write command name starts with "Par" 
  // (i.e. parallel output cases)
  if (StringOps::startsWith(m_writeSolutionStr, "Par")) {
    m_appendRank = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGWriter::openImpl()
{
  CFAUTOTRACE;
  computeFullOutputName();
  m_data->setFilename(m_fullOutputName);
}

//////////////////////////////////////////////////////////////////////////////

void CGWriter::writeImpl()
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

void CGWriter::closeImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::string CGWriter::getFormatExtension() const
{
  CFAUTOTRACE;
  return std::string(".cgns");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD
