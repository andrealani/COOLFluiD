#include <iostream>
#include <fstream>

#include "Environment/CFEnv.hh"
#include "GETModel/GET_DM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::GETModel;

//////////////////////////////////////////////////////////////////////////////

class AppOptions : public Config::ConfigObject {

  public:

    static void defineConfigOptions(Config::OptionList& options)
    {
      options.addConfigOption< CFuint > ("log","Logging level for output");
      options.addConfigOption< bool >   ("help", "Show this help");
    }

    AppOptions() : ConfigObject("Test for Analytical Model lib"),
      logLevel(INFO),
      showUsage(false)
    {
      addOptionsTo(this);
      addEnvironmentVariable("COOLFLUID_OPTS");

      const bool notNested = false;

      setParameter("log",    &logLevel,    notNested);
      setParameter("help",   &showUsage,   notNested);
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
  Common::CFLogger::getInstance().setMainLoggerLevel(options.logLevel);

  // showUsage and exit if any argument is unprocessed
  unusedArgs = options.getOptionList().getUnprocessedArgs();
  if (unusedArgs.size()!=0)
  {
    ConfigArgs::iterator itr = unusedArgs.begin();
    for(; itr != unusedArgs.end(); ++itr) {
      clog << "Unused argument: " << itr->first << " " << itr->second << endl;
    }
    cout << options.writeUsage();
    exit(1);
  }

  // showUsage if --help
  if (options.showUsage)
  {
    cout <<  options.writeUsage();
    exit(0);
  }

  int return_value = 0;
  try
  {
    DomainModel * dm = new GETModel::GET_DM("GET_DM");
    dm->setConfigFile("config.conf");
    dm->configure();

	cout << "after configure" << endl;

    DomainModel::TRidx idx = 1;
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

    pcoord [XX] =  0.75;   // initial guess
    xcoord [XX] =  0.25;
    xcoord [YY] =  0.0594124;
    // 0.497243 -0.0531107
	//0.25 0.0594124

    cout << "XY [" << xcoord << "] : UV initial [" << pcoord << "]\n" << endl;
    dm->computeParamCoord(idx,xcoord,pcoord);
    cout << "XY [" << xcoord << "] : UV calc    [" << pcoord << "]\n" << endl;
    print_data(dm,idx,pcoord,xcoord,deriv);

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
