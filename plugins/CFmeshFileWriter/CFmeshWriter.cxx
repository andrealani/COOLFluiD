// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFmeshFileWriter/CFmeshWriter.hh"
#include "Environment/ObjectProvider.hh"
#include "CFmeshFileWriter/CFmeshFileWriter.hh"
#include "Environment/DirPaths.hh"
#include "Framework/PathAppender.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CFmeshWriter,
               OutputFormatter,
               CFmeshFileWriterModule,
               1>
cfmeshWriterProvider("CFmesh");

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("WriteSol","Command to wrtie the solution.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

CFmeshWriter::CFmeshWriter(const std::string& name)
  : OutputFormatter(name)
{
   addConfigOptionsTo(this);
  _data.reset(new CFmeshWriterData(this));

  _setupStr = "StdSetup";
   setParameter("SetupCom",&_setupStr);

  _unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&_unSetupStr);

   _writeSolutionStr = "ParWriteSolution";
   setParameter("WriteSol",&_writeSolutionStr);
}

//////////////////////////////////////////////////////////////////////////////

CFmeshWriter::~CFmeshWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> CFmeshWriter::getMethodData () const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::setMethodImpl()
{
  // setup parent class
  OutputFormatter::setMethodImpl();

  setupCommandsAndStrategies();
  cf_assert(_setup.isNotNull());
  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::unsetMethodImpl()
{
  cf_assert(_unSetup.isNotNull());
  _unSetup->execute();
  unsetupCommandsAndStrategies();

  // unsetup parent class
  OutputFormatter::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::configure ( Config::ConfigArgs& args )
{
  OutputFormatter::configure(args);
  
  _data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( _data.getPtr(), args );
  
  // add configures to the CFmeshWriterCom's

  configureCommand<CFmeshWriterData,CFmeshWriterComProvider>(args,_setup,_setupStr,_data);

  configureCommand<CFmeshWriterData,CFmeshWriterComProvider>(args,_unSetup,_unSetupStr,_data);

  configureCommand<CFmeshWriterData,CFmeshWriterComProvider>(args,_writeSolution,_writeSolutionStr,_data);
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::openImpl()
{
  computeFullOutputName();
  _data->setFilename(m_fullOutputName);
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::writeImpl()
{
  cf_assert(_writeSolution.isNotNull());
  _writeSolution->execute();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::closeImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string CFmeshWriter::getFormatExtension() const
{
  return std::string(".CFmesh");
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriter::computeFullOutputName()
{
  CFAUTOTRACE;
  
  using namespace boost::filesystem;
  
  path fpath = Environment::DirPaths::getInstance().getResultsDir() / m_filename;
  fpath = PathAppender::getInstance().appendAllInfo(fpath, m_appendIter, m_appendTime, false);
  m_fullOutputName = boost::filesystem::change_extension(fpath, getFormatExtension());
  
  SimulationStatus::getInstance().setLastOutputFile
    (m_fullOutputName,getNamespace(),SubSystemStatusStack::getActive()->getSubSystemName());
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD
