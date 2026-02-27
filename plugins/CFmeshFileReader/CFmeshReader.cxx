// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Common/PEFunctions.hh"
#include "Common/Stopwatch.hh"
#include "Common/CFLog.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/DirPaths.hh"

#include "Framework/SubSystemStatus.hh"
#include "CFmeshFileReader/CFmeshReader.hh"
#include "CFmeshFileReader/CFmeshFileReader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

ObjectProvider<CFmeshReader,
               MeshCreator,
               CFmeshFileReaderPlugin,
               1>
CFmeshReaderProvider("CFmeshFileReader");

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string > ("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string > ("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string > ("ReadCFmesh","Command to read the CFmesh. This command seldomly needs overriding.");
   options.addConfigOption< std::string > ("convertFrom","Name of format from which to convert to CFmesh.");
   options.addConfigOption< bool > ("convertBack","Also convert back to the original format. Usefull only for debugging.");
   options.addConfigOption< bool > ("onlyConversion","Only convert the mesh without loading it into memory.");
}

//////////////////////////////////////////////////////////////////////////////

CFmeshReader::CFmeshReader(const std::string& name)  : MeshCreator(name)
{
  addConfigOptionsTo(this);
  m_data.reset(new CFmeshReaderData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&m_unSetupStr);

/// @todo this should be removed when ParCFmeshFile reader can
///       work also in serial mode
#ifdef CF_HAVE_MPI
  m_readCFmeshStr = "ParReadCFmesh";
#else
  m_readCFmeshStr = "StdReadCFmesh";
#endif // CF_HAVE_MPI
  setParameter("ReadCFmesh",&m_readCFmeshStr);

  m_converterStr = "Null";
  setParameter("convertFrom",&m_converterStr);

  m_onlyConversion = false;
  setParameter("onlyConversion",&m_onlyConversion);

  m_convertBack = false;
  setParameter("convertBack",&m_convertBack);
}

//////////////////////////////////////////////////////////////////////////////

CFmeshReader::~CFmeshReader()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> CFmeshReader::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::generateMeshDataImpl()
{
  CFAUTOTRACE;

  // convert only if m_converterStr != "CFmesh", which assumes
  // that you create the mesh from CFmesh format
  
  CFLog(VERBOSE, "CFmeshReader::generateMeshDataImpl() => start\n"); 
  
  if (m_converterStr != "CFmesh")
  {
    convertFormat();
    CFLog(VERBOSE,"CFmeshReader : finished converting\n");
  }
 
  if (!m_onlyConversion) {
    Common::Stopwatch<WallTime> stp;
    CFLog(VERBOSE, "CFmeshReader::generateMeshDataImpl() => 2\n"); 
    
    stp.start();
    
    cf_assert(m_readCFmesh.isNotNull());
    m_readCFmesh->execute();
    
    stp.stop();
    
    CFLog(NOTICE, "Building the mesh took: " << stp.read() << "s\n");
  }
  
  CFLog(VERBOSE, "CFmeshReader::generateMeshDataImpl() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::processMeshDataImpl()
{
  CFAUTOTRACE;
  //There is nothing to be done
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::buildMeshDataImpl()
{
  CFAUTOTRACE;
  //There is nothing to be done
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::convertFormat()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "CFmeshReader::convertFormat() => start\n");
  
  // build the MeshFormatConverter to perform the conversion
  
  // AL: this trick is necessary because otherwise the corresponding Converter
  // options will be removed from the args after this delayed configuration
  // originalArgs != m_stored_args after a call to configureNested() 
  // To be investigated whether this is a bug...
  Config::ConfigArgs originalArgs = m_stored_args;

  Common::SelfRegistPtr<MeshFormatConverter>* converter =
   FACTORY_GET_PROVIDER(getFactoryRegistry(), MeshFormatConverter, m_converterStr)->createPtr(m_converterStr);
 
  CFLog(VERBOSE, "CFmeshReader::convertFormat() => 2\n");

  cf_assert(converter->isNotNull());
  (*converter)->setFactoryRegistry(getFactoryRegistry());
  configureNested ( converter->getPtr(), originalArgs );

  CFLog(VERBOSE, "CFmeshReader::convertFormat() => 3\n");

  const std::string nsp = getMethodData()->getNamespace();
  
  // ensure the converter works serially only on the processor with rank 0
  // if it is not enabled to work in parallel
  if (PE::GetPE().IsParallel() && !(*converter)->isParallel()) {
    PE::GetPE().setBarrier(nsp);
    if (PE::GetPE().GetRank(nsp) == 0) { convert(*converter); }
    PE::GetPE().setBarrier(nsp);
  }
  else { 
    convert(*converter);
  }
   
  delete converter;
 
  CFLog(VERBOSE, "CFmeshReader::convertFormat() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
void CFmeshReader::convert(Common::SelfRegistPtr<MeshFormatConverter> converter)
{
  CFAUTOTRACE;
  
  using namespace boost::filesystem;

  path fromfile = DirPaths::getInstance().getWorkingDir() / m_data->getConvertFromFileName();
  path tofile   = DirPaths::getInstance().getWorkingDir() / m_data->getFileName();

#ifndef NDEBUG
  converter->checkFormat(fromfile);
#endif

  converter->convert(fromfile, tofile);

  if (m_convertBack)
  {
    path backfile(DirPaths::getInstance().getWorkingDir());
    backfile /= m_data->getFileName();
#ifdef CF_HAVE_BOOST_1_85
    backfile.replace_extension("-convertedback" + converter->getOriginExtension());
#else
    backfile = boost::filesystem::change_extension(backfile, "-convertedback" + converter->getOriginExtension());
#endif  
    converter->convertBack(backfile);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::modifyFileNameForRestart(const std::string filename)
{
  m_data->setFileName(filename);
  m_data->setIsTranslated(false);
  m_converterStr = "CFmesh";
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::modifyFileName(const std::string filename)
{
  m_data->setFileName(filename);
  m_data->setConvertFromFileName(filename);
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  MeshCreator::configure(args);
  
  m_data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_data.getPtr(), args );
  
  configureCommand<CFmeshReaderData,
                   CFmeshReaderComProvider>(args,
                                            m_setup,
                                            m_setupStr,
                                            m_data);

  configureCommand<CFmeshReaderData,
                   CFmeshReaderComProvider>(args,
                                            m_unSetup,
                                            m_unSetupStr,
                                            m_data);

  configureCommand<CFmeshReaderData,
                   CFmeshReaderComProvider>(args,
                                            m_readCFmesh,
                                            m_readCFmeshStr,
                                            m_data);

  m_stored_args = args;
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::setMethodImpl()
{
  setupCommandsAndStrategies();
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReader::unsetMethodImpl()
{
  cf_assert(m_unSetup.isNotNull());
  m_unSetup->execute();
  unsetupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
