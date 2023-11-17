#include "FiniteVolume/DistanceBasedExtrapolatorMagnetogram.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/OSystem.hh"
#include "Common/StringOps.hh" 

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolatorMagnetogram,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeModule>
distanceBasedExtrapolatorMagnetogramProvider("DistanceBasedMagnetogram");

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorMagnetogram::DistanceBasedExtrapolatorMagnetogram
(const std::string& name) :
  DistanceBasedExtrapolator<CellCenterFVMData>(name)
{
  addConfigOptionsTo(this);

  _sigma = 8.; // between 8 and 10
  setParameter("Sigma",&_sigma);

  _scaling_factor = 0.26;
  setParameter("ScalingFactor",&_scaling_factor);
  
  _link = "https://gong.nso.edu/data/magmap/QR/bqj/201207/mrbqj120712/mrbqj120712t1154c2125_024.fits.gz";
  setParameter("Link",&_link);
  
  _pyCommand = "python3";
  setParameter("PyCommand",&_pyCommand);
  
  _Brefval = 2.2e-4;
  setParameter("Brefval",&_Brefval);
  
  _runPyScript = true;
  setParameter("RunPyScript",&_runPyScript);
}
      
//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorMagnetogram::~DistanceBasedExtrapolatorMagnetogram()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorMagnetogram::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Sigma","Smoothing width in theta/phi");
  options.addConfigOption< CFreal >("ScalingFactor","Scale up/down the magnetogram B-field strength");
  options.addConfigOption< string >("Link","The URL for downloading a specific magnetogram");
  options.addConfigOption< string >("PyCommand", "Python command to run to prepare the magnetogram.");
  options.addConfigOption< CFreal >("Brefval","Reference value (conversion to adimensional units for the B-field");
  options.addConfigOption< bool >("RunPyScript","Flag telling whether to run the Python script to prepare the magnetogram.");
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorMagnetogram::setup()
{
  static bool firstIn = true;
  cf_assert(m_fileNameTw.size() == 1);
  if (firstIn) {
    if (_runPyScript) { 
      CFLog(INFO, "DistanceBasedExtrapolatorMagnetogram::setup() => PREPARE MAGNETOGRAM\n");
      
      // Run a python script for preparing the magnetogram:
      std::string cmd = _pyCommand + " prepare-the-magnetogram.py " + _link + " " +
	StringOps::to_str(_sigma) + " " + StringOps::to_str(_scaling_factor) +  " " + 
	StringOps::to_str(_Brefval) + " ; cp -r magnetogram*.dat " + m_fileNameTw[0];
      
      CFLog(INFO, "DistanceBasedExtrapolatorMagnetogram::setup() => running \"" << cmd << "\"\n");
      
      Common::OSystem::getInstance().executeCommand(cmd);
      firstIn = false;
    }
  }
  
  DistanceBasedExtrapolator<CellCenterFVMData>::setup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
