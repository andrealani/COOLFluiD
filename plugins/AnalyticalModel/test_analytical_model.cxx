// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>
#include <fstream>



#include "Environment/CFEnv.hh"
#include "AnalyticalModel/AnalyticalDM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::AnalyticalModel;

//////////////////////////////////////////////////////////////////////////////

class AppOptions : public Config::ConfigObject {

  public:

    static void defineConfigOptions(Config::OptionList& options)
    {
      options.addConfigOption< CFuint > ("log","Logging level for output");
      options.addConfigOption< bool >   ("help", "Show this help");
    }

    AppOptions() : ConfigObject("TestAnalyticalModel"),
      logLevel(INFO),
      showUsage(false)
    {
      addConfigOptionsTo(this);
      addEnvironmentVariable("COOLFLUID_OPTS");

      setParameter("log",    &logLevel);
      setParameter("help",   &showUsage);
    }

  public:

      CFuint  logLevel;
      bool showUsage;
};

//////////////////////////////////////////////////////////////////////////////

void print_data(DomainModel * dm,
                const DomainModel::TRidx& idx,
                const RealVector& pcoord,
                RealVector& xcoord,
                std::vector<RealVector>& deriv)
{
  dm->computeCoord(idx,pcoord,xcoord);
  dm->compute1stDeriv(idx,pcoord,deriv);
  cout << "UV [" << pcoord << "] : XY [" << xcoord << "] : DXYDUV [" << deriv[0] << "] \n" << endl;
}

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  AppOptions options;
  options.getOptionList().setStrictArgs(true);

  // unused arguments
  ConfigArgs unusedArgs;

  // process the command line
  try
  {
    options.getOptionList().processCommandLine(argc, argv);
  }
  // Something bad happened. Dump the usage statement:
  catch (const Exception& e)
  {
    cerr << options.writeUsage();
    exit(1);
  }

  Environment::CFEnv::getInstance().initiate(argc, argv);
  CFLogger::getInstance().setMainLoggerLevel(options.logLevel);

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

  // showUsage if --help
  if (options.showUsage)
  {
    cout << options.writeUsage();
    exit(0);
  }

//   Environment::DirPaths::getInstance().setBaseDir(options.baseDir);
//   Environment::DirPaths::getInstance().addModuleDirs(options.libDir);

  int return_value = 0;
  try
  {
    DomainModel * dm = new AnalyticalModel::AnalyticalDM("AnalyticalDM");
    dm->setConfigFile("config.conf");
    ConfigArgs args;
    dm->configure(args);

    DomainModel::TRidx idx = 0;
    RealVector pcoord (DIM_1D);
    RealVector xcoord (DIM_2D);
    std::vector<RealVector> deriv (DIM_1D,xcoord);

    pcoord[XX] = 0.0;
    print_data(dm,idx,pcoord,xcoord,deriv);

    pcoord[XX] = 0.25;
    print_data(dm,idx,pcoord,xcoord,deriv);

    pcoord[XX] = 0.5;
    print_data(dm,idx,pcoord,xcoord,deriv);

    pcoord[XX] = 0.75;
    print_data(dm,idx,pcoord,xcoord,deriv);

    pcoord[XX] = 1.0;
    print_data(dm,idx,pcoord,xcoord,deriv);

    pcoord [XX] =  0.75;
    xcoord [XX] =  7.07107;
    xcoord [YY] =  7.07107;

    cout << "XY [" << xcoord << "] : UV [" << pcoord << "]\n" << endl;
    dm->computeParamCoord(idx,xcoord,pcoord);
    cout << "XY [" << xcoord << "] : UV [" << pcoord << "]\n" << endl;

  }
  catch (std::exception& e)
  {
    cerr << e.what() << endl;
    cerr << "Aborting ... " << endl;
    return_value = 1;
  }
  catch (...)
  {
    cerr << "Unknown exception thrown and not caught !!!" << endl;
    cerr << "Aborting ... " << endl;
    return_value = 1;
  }

  Environment::CFEnv::getInstance().terminate();

  return return_value;
}
