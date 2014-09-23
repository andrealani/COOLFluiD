// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>
#include <iterator>

#include <boost/filesystem/path.hpp>

#include "Common/xmlParser.h"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Config/ConfigObject.hh"
#include "Config/ConfigFileReader.hh"
#include "Config/ConverterTools.hh"
#include "Config/XMLConfigFileReader.hh"

#include "Environment/CFEnv.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

// template<typename T>
// void print_vector (std::ostream& ostr, const std::vector<T>& t, const std::string& delimiter)
// {
//   std::copy(t.begin(), t.end(), std::ostream_iterator<T>(ostr, delimiter));
// }
//
// template < typename T, typename InputIterator >
// void print (std::ostream& ostr,
//             InputIterator itbegin,
//             InputIterator itend,
//             const std::string& delimiter)
// {
//   std::copy(itbegin, itend, std::ostream_iterator<T>(ostr, delimiter));
// }

//////////////////////////////////////////////////////////////////////////////

XMLNode getNode ( XMLNode base_node, const std::vector<std::string>& hierarchy)
{
  XMLNode child = base_node;
  for (std::vector<std::string>::const_iterator it = hierarchy.begin(),
       stop = hierarchy.end();
       it != stop; ++it)
  {
    const std::string name (*it);
    XMLNode xNode = child.getChildNode( name.c_str(), 0 ); // assuming only no duplicated names in each node
    if (xNode.isEmpty()) // create this node
    {
      xNode = child.addChild(name.c_str());
    }
    child = xNode;
  }
  return child;
}

//////////////////////////////////////////////////////////////////////////////

class AppOptions : public Config::ConfigObject {

  public:

    static void defineConfigOptions(Config::OptionList& options)
    {
      options.addConfigOption< std::string > ("cfcase","Name of CFcase file");
      options.addConfigOption< std::string > ("xcfcase","Name of XCFcase file");
      options.addConfigOption< bool >     ("backward", "Convert from XCFcase to CFcase");
      options.addConfigOption< CFuint >   ("log","Logging level for output");
      options.addConfigOption< std::string > ("bdir", "The base dir to preced all file and dir paths.");
      options.addConfigOption< bool >     ("help", "Show this help");
      options.addConfigOption< std::vector<std::string> >("ldir", "The list of dir paths where to search for module plug-in libraries. Paths are fully qualified and not preceeded with base dir.");
    }

    AppOptions() : ConfigObject("XCFcaseConverter"),
      baseDir(""),
      libDir(0),
      cfcase(""),
      xcfcase(""),
      backward(false),
      logLevel(INFO),
      showUsage(false)
    {
      addConfigOptionsTo(this);
      addEnvironmentVariable("XCFCASE_CONVERTER_OPTS");

      setParameter("cfcase",   &cfcase);
      setParameter("xcfcase",  &xcfcase);
      setParameter("backward", &backward);
      setParameter("log",      &logLevel);
      setParameter("bdir",     &baseDir);
      setParameter("ldir",     &libDir);
      setParameter("help",     &showUsage);
    }

  public:

      std::string baseDir;
      std::vector<std::string> libDir;
      std::string cfcase;
      std::string xcfcase;
      bool backward;
      CFuint  logLevel;
      bool showUsage;
};

//////////////////////////////////////////////////////////////////////////////

struct LabelValue
{
  LabelValue ( const std::pair< std::string, std::string >& arg )
   : value(arg.second)
  {
    nested_labels = StringOps::getWords(arg.first,'.');
    label = nested_labels.back();
    nested_labels.pop_back();

//     cout << "nested_labels [" ;
//     std::copy(nested_labels.begin(), nested_labels.end(), std::ostream_iterator<std::string>(std::cout, ","));
//     cout << "]    \tlabel [" << label ;
//     cout << "]    \tvalue [" << value << "]" << std::endl;
  }

  std::vector<std::string> nested_labels;
  std::string              label;
  std::string              value;
};

//////////////////////////////////////////////////////////////////////////////

ConfigArgs * load_cfcase(const std::string& name)
{
    // to be deleted
    ConfigFileReader configFile;
    ConfigArgs * args = new ConfigArgs();
    configFile.parse(name, *args);

    cout << "loading file [" << name << "]" << endl;

//     cout << args.str() << endl;
    return args;
}

//////////////////////////////////////////////////////////////////////////////

void write_cfcase(const std::string& cfcase, const ConfigArgs& args)
{
    Common::SelfRegistPtr<Environment::FileHandlerOutput> filehandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    cf_assert (!filehandle->isopen());
    boost::filesystem::path filepath = cfcase;
    ofstream& fout = filehandle->open(filepath);

    cout << "saving file [" << cfcase << "]" << endl;

    for (ConfigArgs::const_iterator it = args.begin(), stop = args.end(); it != stop; ++it)
    {
      fout << it->first << " = " << it->second << "\n";
    }
}

//////////////////////////////////////////////////////////////////////////////

void write_xcfcase(const std::string& xcfcase, const ConfigArgs& args)
{
 
 Common::SelfRegistPtr<Environment::FileHandlerOutput> filehandle =
   Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
 cf_assert (!filehandle->isopen());
 boost::filesystem::path filepath = xcfcase;
 ofstream& fout = filehandle->open(filepath);
 
 XMLNode xMainNode = ConverterTools::configArgsToXCFcase(args);

 cout << "saving file [" << xcfcase << "]" << endl;

    // create a string to hold the file contents
 char* tstr = xMainNode.createXMLString();
 fout << tstr << endl;
 free(tstr);
}

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  bool return_value = true;
  AppOptions options;
  options.getOptionList().setStrictArgs(true);

  // unused arguments
  ConfigArgs unusedArgs;

  // process the command line
  try { options.getOptionList().processCommandLine(argc, argv); }
  // Something bad happened. Dump the usage statement:
  catch (const Exception& e)
  {
    cerr << options.writeUsage();
    exit(1);
  }

  // set the command line log level
  // takes precedence over the config
  CFLogger::getInstance().setMainLoggerLevel(options.logLevel);

  // build the environment
  CFEnv& cf_env = CFEnv::getInstance();
  // initiate the environemnt
  cf_env.initiate(argc, argv);

  // set essential paths
  DirPaths::getInstance().setBaseDir(options.baseDir);
  DirPaths::getInstance().addModuleDirs(options.libDir);

  // configure the environment
  ConfigArgs env_args;
  cf_env.configure(env_args);

  // setup the runtime environment
  cf_env.setup();

  // showUsage and exit if any argument is unprocessed
  /// @todo when processing command line detect which options where not used
//   unusedArgs = options.getOptionList().getUnprocessedArgs();
//   if (unusedArgs.size()!=0)
//   {
//     ConfigArgs::iterator itr = unusedArgs.begin();
//     for(; itr != unusedArgs.end(); ++itr) {
//       clog << "Unused argument: " << itr->first << " " << itr->second << endl;
//     }
//     cout << options.writeUsage();
//     exit(1);
//   }

  // no cfcase selected
  if(options.cfcase == "")
  {
    options.showUsage = true;
    clog << "No origin CFcase file chosen!" << endl;
  }

  // no xcfcase selected
  if(options.xcfcase == "")
  {
    options.showUsage = true;
    clog << "No destination XCFcase file chosen!" << endl;
  }

  // showUsage if --help
  if (options.showUsage)
  {
    cout << options.writeUsage();
    exit(0);
  }

  // real stuff begins
  try
  {
    ConfigArgs * args;

    if (!options.backward)
    {
      args = load_cfcase(options.cfcase);
      // write_cfcase("output.CFcase",*args);
      write_xcfcase(options.xcfcase,*args);
    }
    else
    {
      args = load_cfcase(options.xcfcase);
      write_cfcase(options.cfcase,*args);
      // write_xcfcase("outputx.XCFcase",*args);
    }

      deletePtr(args);
  }
  catch (std::exception& e) {
    cerr << e.what() << endl;
    cerr << "Aborting ... " << endl;
    return_value = 1;
  }
  catch (...) {
    cerr << "Unknown exception thrown and not caught !!!" << endl;
    cerr << "Aborting ... " << endl;
    return_value = 1;
  }

  // unsetup the runtime environment
  cf_env.unsetup();
  // terminate the runtime environment
  cf_env.terminate();

  return return_value;
}
