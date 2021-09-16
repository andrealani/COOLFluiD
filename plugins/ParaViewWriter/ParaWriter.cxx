// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ParaWriter.hh"
#include "Environment/ObjectProvider.hh"
#include "ParaViewWriter/ParaViewWriter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ParaWriter,
               OutputFormatter,
               ParaViewWriterModule,
               1>
ParaWriterProvider("ParaView");

//////////////////////////////////////////////////////////////////////////////

void ParaWriter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("WriteSol","Command to write the solution.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

ParaWriter::ParaWriter(const std::string& name)
  : OutputFormatter(name)
{
  addConfigOptionsTo(this);
  m_data.reset(new ParaWriterData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

  m_writeSolutionStr = "WriteSolutionNoOverlap"; // "WriteSolution";
  setParameter("WriteSol",&m_writeSolutionStr);
}

//////////////////////////////////////////////////////////////////////////////

ParaWriter::~ParaWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> ParaWriter::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void ParaWriter::setMethodImpl()
{
  OutputFormatter::setMethodImpl();

  setupCommandsAndStrategies();
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void ParaWriter::unsetMethodImpl()
{
  cf_assert(m_unSetup.isNotNull());
  m_unSetup->execute();
  unsetupCommandsAndStrategies();

  OutputFormatter::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void ParaWriter::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  OutputFormatter::configure(args);

  m_data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_data.getPtr(), args );

  // add configures to the ParaWriterCom's

  configureCommand<ParaWriterData,ParaWriterComProvider>( args, m_setup, m_setupStr, m_data);

  configureCommand<ParaWriterData,ParaWriterComProvider>( args, m_unSetup, m_unSetupStr, m_data);

  configureCommand<ParaWriterData,ParaWriterComProvider>( args, m_writeSolution, m_writeSolutionStr, m_data);
}

//////////////////////////////////////////////////////////////////////////////

void ParaWriter::openImpl()
{
  CFAUTOTRACE;
  computeFullOutputName();
  m_data->setFilename(m_fullOutputName);
}

//////////////////////////////////////////////////////////////////////////////

void ParaWriter::writeImpl()
{
  CFAUTOTRACE;
  cf_assert(isSetup());
  cf_assert(isConfigured());

  cf_assert(m_writeSolution.getPtr() != CFNULL);
  // set the nodal states at first
  m_data->getCollaborator<SpaceMethod>()->extrapolateStatesToNodes();
  // write the solution file
  m_writeSolution->execute();
}

//////////////////////////////////////////////////////////////////////////////

void ParaWriter::closeImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::string ParaWriter::getFormatExtension() const
{
  CFAUTOTRACE;
  return std::string(".vtu");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD
