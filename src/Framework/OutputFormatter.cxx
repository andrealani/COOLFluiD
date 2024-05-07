// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/filesystem/convenience.hpp>

#include "Framework/OutputFormatter.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/PathAppender.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("FileName","The name with extension of the file to write the output to.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >("SaveRate","Rate to intermediate solution to Tecplot file.");
  options.addConfigOption< bool >("SaveFinal", "Save again on finishing simulation.");
  options.addConfigOption< bool >("AppendTime","Save each iteration to different file with suffix m_time#.");
  options.addConfigOption< bool >("AppendIter","Save each iteration to different file with suffix m_iter#.");
  options.addConfigOption< bool >("AppendRank","Append the processor rank to the file.");
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::build_dynamic_functions()
{
  add_dynamic_function("open",&OutputFormatter::open);
  add_dynamic_function("close",&OutputFormatter::close);
  add_dynamic_function("write",&OutputFormatter::write);
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::registActionListeners()
{
  Method::registActionListeners();

  // const std::string ssname = SubSystemStatusStack::getCurrentName();
  // event_handler->addListener(event_handler->key(ssname, "CF_ON_OUTPUTFORMATTER_WRITE"),
  //                            this,&OutputFormatter::write);
}

//////////////////////////////////////////////////////////////////////////////

OutputFormatter::OutputFormatter(const std::string& name) : Method(name)
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);

  setParameter("FileName",&m_filenameStr);

  m_saveRate = 100;
  setParameter("SaveRate",&m_saveRate);

  m_savefinal = true;
  setParameter("SaveFinal",&m_savefinal);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
  
  m_appendRank = true;
  setParameter("AppendRank",&m_appendRank);
}

//////////////////////////////////////////////////////////////////////////////

OutputFormatter::~OutputFormatter()
{
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::setMethodImpl()
{
  if ((m_filenameStr.empty()) && (!isNonRootMethod()))
  {
    throw Common::FilesystemException (FromHere(),"OutputFormatter: no output file specified");
  }

  // remove the extension and assign to filename
  m_filename = boost::filesystem::basename(boost::filesystem::path(m_filenameStr));

  SubSystemStatusStack::getActive()->setAppendToFile(m_appendIter,m_appendTime);
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::setOutputFileName(const std::string filename)
{
  // remove the extension and assign to filename
  m_filename = boost::filesystem::basename(boost::filesystem::path(filename));
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::computeFullOutputName()
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
  
  path fpath = Environment::DirPaths::getInstance().getResultsDir() / m_filename;
  fpath = PathAppender::getInstance().appendAllInfo(fpath, m_appendIter,m_appendTime,m_appendRank);
  m_fullOutputName = boost::filesystem::change_extension(fpath, getFormatExtension());
}
    
//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::open()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  openImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::write()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  writeImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void OutputFormatter::close()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  closeImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

bool OutputFormatter::isSaveNow( const bool force_write )
{
  // force writing happens iff save final is OK
  if ( force_write ) return m_savefinal;
  
  // no save rate means never save must return before divide by zero
  if (m_saveRate == 0) return false;
  
  NamespaceSwitcher& nsw = NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName());
  const string sssName = nsw.getName(mem_fun<string,Namespace>(&Namespace::getSubSystemStatusName), true);
  SafePtr<SubSystemStatus> currSSS = SubSystemStatusStack::getInstance().getEntry(sssName);
  cf_assert(currSSS.isNotNull());
  
  // if divide by save rate is integer then it is time to save
  if (!(currSSS->getNbIter() % m_saveRate)) {
    return true;
  }
  else {
    return false;
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
