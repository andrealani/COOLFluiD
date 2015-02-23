// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdlib> // for getenv
#include <fstream>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>

#include "Common/CFLog.hh"
#include "Common/PE.hh"
#include "Common/FilesystemException.hh"

#include "Config/ConfigObject.hh"
#include "Config/ConfigFileReader.hh"
#include "Config/ConfigFacility.hh"
#include "Config/XMLConfigFileReader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Config {

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct indent : public std::unary_function<T, void>
{
  indent(std::ostream& out) : os(out) {}
  void operator()(T x) { os << "    " << x << "\n"; }
  std::ostream& os;
};

//////////////////////////////////////////////////////////////////////////////

const std::string ConfigObject::NEST_SEPARATOR = ".";
const std::string ConfigObject::NEST_SEPARATOR_ENVVAR = "_";

//////////////////////////////////////////////////////////////////////////////

ConfigObject::ConfigObject(const std::string& name,
                           const bool&   helpOnError,
                           const std::string& usageSummary)
  : NamedObject(name),
    m_options(),
    m_nested_objs(),
    m_config(false),
    m_helpOnError(helpOnError),
    m_hadError(false),
    m_registered(false),
    m_config_file(),
    m_env_vars(),
    m_nest(),
    m_usageSummary(usageSummary),
    parentConfig ( CFNULL )
{
  cf_assert(name != "");
}

//////////////////////////////////////////////////////////////////////////////

ConfigObject::~ConfigObject()
{
 unregistFromConfigFacility();

 // removing object from parent vector
 if(parentConfig != CFNULL)
 {
  parentConfig->unregisterFromParent(this);
  parentConfig = CFNULL;
 }
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::configure_xml ( const std::string& xml_args )
{
	// convert from xml to ConfigArgs
    ConfigArgs config_args;
    XMLConfigFileReader xmlConfig;
    xmlConfig.parseString(xml_args, config_args);

	// call the canonical configure() function
	this->configure(config_args);
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::configure ( Config::ConfigArgs& args )
{
  m_config = true;

  std::string basename = getNestName() + ConfigObject::NEST_SEPARATOR;
  // prepend to each option the basename so they can be matched
  m_options.setHierarchyOptions ( basename );

  // process the config arguments
  processOptions ( args );

#ifndef NDEBUG
  std::string debug_info = m_options.debugInfo();

  CFuint rank = 0;
  if ( Common::PE::IsInitialised() )
    rank = Common::PE::GetPE().GetRank();
  if ( rank == 0 )
  {
    ofstream configdbg ( "config-debug-info.log" ,ios_base::app);
    configdbg << "### " << getNestName() << "\n";
    configdbg << debug_info << "\n";
  }
#endif

  // remove the base name form all options
  m_options.setHierarchyOptions ( std::string() );
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::configureDynamic ( ConfigArgs& args )
{
  std::string basename = getNestName() + ConfigObject::NEST_SEPARATOR;
  // prepend to each option the basename so they can be matched
  m_options.setHierarchyOptions ( basename );

  // do dynamic configuration
  bool dynamic = true;
  m_options.processArgs ( args, dynamic );

  // remove the base name form all options
  m_options.setHierarchyOptions ( std::string() );
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::configureNested ( ConfigObject* nestConfig, ConfigArgs& args )
{
  cf_assert(nestConfig != CFNULL);
  configureNested(*nestConfig, args);
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::configureNested (ConfigObject& nestConfig, ConfigArgs& args )
{
  addNestedConfigObject(nestConfig);

  nestConfig.setNest ( getNestName() );      // sets the name of the nesting object in this object
  nestConfig.configure( args );  // do the configuration configure
}

//////////////////////////////////////////////////////////////////////////////

Option * ConfigObject::getOption(const std::string& name) const
{
  // gets the option
  Option * opt = m_options.findFirstByName(name);
  if (opt == CFNULL)
    throw ConfigOptionException (FromHere(), "No option exists with name : " + name);
  return opt;
}


//////////////////////////////////////////////////////////////////////////////

void ConfigObject::unconfigure()
{
  m_config = false;
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::setNest (const std::string& nest)
{
  m_nest = nest;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConfigObject::getNest() const
{
  return m_nest;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConfigObject::getNestName() const
{
  return m_nest.empty() ? getName() : m_nest +  ConfigObject::NEST_SEPARATOR + getName();
}

//////////////////////////////////////////////////////////////////////////////

std::string ConfigObject::usageSummary() const
{
  return m_usageSummary;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConfigObject::writeUsage() const
{
  std::ostringstream os;

  os << "Usage: " << getName() << " " << usageSummary() << "\n";
  os << "\n";
  os << "Options:\n";
  const OptionList& options = getOptionList();
  os << options.writeUsage();

  const std::vector<std::string>& env_vars = getEnvVars();
  if (!env_vars.empty())
  {
    os << "\nEnvironment Variables\n";
    for_each(env_vars.begin(), env_vars.end(), indent< std::string >(os));
  }

  std::string conffile ( getConfigFile().string() );
  if (!conffile.empty())
  {
    os << "\nConfiguration File:" << conffile  << "\n";
  }
  return os.str();
}

//////////////////////////////////////////////////////////////////////////////

std::string ConfigObject::debugInfo () const
{
  std::ostringstream oss;

  const OptionList& options = getOptionList();

//   ConfigArgs unprArgs = options.getUnprocessedArgs();
//   ConfigArgs procArgs = options.getProcessedArgs();

  oss << "#######################################################\n";
  oss << "# Debug Config Info: " << getName() << "\n";
  oss << "#  OptionList:\n";
  oss << options.debugInfo() << "\n";
  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::processConfigFile ( ConfigArgs& args )
{
  if( !m_config_file.string().empty() ) // only parse if filename was set
  {

    if ( !boost::filesystem::exists(m_config_file) || boost::filesystem::is_directory(m_config_file) )
      throw Common::FilesystemException (FromHere(),"Could not open configuration file: " + m_config_file.string());

    ConfigFileReader config_parser;
    config_parser.parse(m_config_file, args);
  }

  m_options.processArgs ( args );
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::processEnvironmentVariables ( ConfigArgs& args )
{
  std::pair<bool, std::string> ret;
  vector<std::string>::iterator it   = m_env_vars.begin();
  vector<std::string>::iterator stop = m_env_vars.end();
  for ( ; it != stop; ++it)
  {
    const std::string& envVar = *it;
    ret = readEnvironmentVariable(envVar);
    // The first member indicates if the second is valid or not
    if (ret.first) { args[envVar] = ret.second; }
  }

  m_options.processArgs ( args );
}

//////////////////////////////////////////////////////////////////////////////

bool ConfigObject::hadError() const
{
  return m_hadError;
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::setConfigFile ( const boost::filesystem::path& file )
{
  m_config_file = file;
}

//////////////////////////////////////////////////////////////////////////////

std::pair<bool,std::string>
ConfigObject::readEnvironmentVariable(const std::string& envVar)
{
  std::string envVarNoDot(envVar);
  Common::StringOps::subst(NEST_SEPARATOR,NEST_SEPARATOR_ENVVAR,envVarNoDot);
  CFchar* ptrValue = getenv(envVarNoDot.c_str());
  if (ptrValue)
  {
    return std::pair<bool,std::string>(true, ptrValue);
  }
  else
  {
    return std::make_pair(false, std::string());
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::addEnvironmentVariable(const std::string& envVar)
{
  m_env_vars.push_back(envVar);
}


//////////////////////////////////////////////////////////////////////////////

void ConfigObject::processParentArgs ( ConfigArgs& args )
{
  std::pair<ConfigArgs,ConfigArgs> procargs = m_options.processArgs ( args );

  // write the config log file
  /// @todo add this in XML CFcase format

  CFuint rank = 0;
  if ( Common::PE::IsInitialised() )
    rank = Common::PE::GetPE().GetRank();

  // only rak 0 writes the configuration file
  if ( rank == 0) {
    std::string filename;
    filename += "config-p" + StringOps::to_str(rank) + ".log";
    ofstream confOut(filename.c_str(),ios_base::app);
    confOut << "### " << getNestName() << "\n";
    confOut << procargs.first.str();
    confOut.close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::processOptions ( ConfigArgs& args )
{
    ConfigArgs& parent_args = args;

    // default does not propagate file config out of the object
    ConfigArgs file_args;
    processConfigFile (file_args);

    // process arguments comming from the parent object
    processParentArgs ( parent_args );

    // default does not propagate envvars config out of the object
    ConfigArgs env_args;
    processEnvironmentVariables (env_args);
    m_options.processArgs ( env_args );

    m_options.doValidation();
    m_options.doComplete();
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::registToConfigFacility()
{
  if ( !m_registered )
    ConfigFacility::getInstance().registConfigObj(this);
  m_registered = true;
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::unregistFromConfigFacility()
{
  if ( m_registered )
    ConfigFacility::getInstance().unregistConfigObj(this);
  m_registered = false;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConfigObject::getTreeXML()
{
//   CFlog << "Dumping XML Tree for ConfigObject: " <<  getName() << " [" <<  getNestName() << "]\n";

  std::string result;

  std::ostringstream os;

  os << "<" << getName() << "\t";
  os << "tree=\"" << "object" << "\" ";
//   os << "type=\"" << DEMANGLED_TYPEID() << "\" ";
  os << "mode=\"basic\" ";
  os << "dynamic=\"static\" ";
  os << ">\n";

  os << m_options.getOptionsXML();

  // loop over nested objects
  std::vector<ConfigObject*>::iterator itr  = m_nested_objs.begin();
  std::vector<ConfigObject*>::iterator last = m_nested_objs.end();
  for (; itr != last; ++itr)
  {
    os << (*itr)->getTreeXML();
  }

  os << "</" << getName() << ">\n";

  return os.str();
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::addNestedConfigObject(ConfigObject& obj)
{
  m_nested_objs.push_back(&obj);

  std::sort(m_nested_objs.begin(), m_nested_objs.end(), std::less<ConfigObject*>());
  std::vector<ConfigObject*>::iterator last_obj = std::unique(m_nested_objs.begin(),m_nested_objs.end());
  m_nested_objs.erase(last_obj, m_nested_objs.end());

  obj.parentConfig = this;
}

//////////////////////////////////////////////////////////////////////////////

void ConfigObject::unregisterFromParent(ConfigObject * config)
{
//   CFlog << "unregistering object [" << config->getName() << "]\n\tfrom parent [" << getName() << "] with [" << m_nested_objs.size() << "] nested objs\n" << CFendl;

  std::sort(m_nested_objs.begin(), m_nested_objs.end(), std::less<ConfigObject*>());
  std::vector<ConfigObject*>::iterator last_obj = std::unique(m_nested_objs.begin(),m_nested_objs.end());
  m_nested_objs.erase(last_obj, m_nested_objs.end());

  std::vector < ConfigObject* >::iterator itr = find(m_nested_objs.begin(), m_nested_objs.end(), config);

  if ( itr != m_nested_objs.end() )
  {
//     CFlog << "removing object [" << config << "] found with iterator [" << *itr << "]\n" << CFendl;
    m_nested_objs.erase(itr);
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Config
} // namespace COOLFluiD
