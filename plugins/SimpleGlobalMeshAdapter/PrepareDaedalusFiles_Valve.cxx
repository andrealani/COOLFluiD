// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "PrepareDaedalusFiles_Valve.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/OSystem.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Common/PE.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PrepareDaedalusFiles_Valve, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> PrepareDaedalusFiles_ValveProvider("PrepareDaedalusFiles_Valve");

//////////////////////////////////////////////////////////////////////////////

void PrepareDaedalusFiles_Valve::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Script","Name of the script to run.");
  options.addConfigOption< std::string >("ScriptStartFile","Name of the startfile for the script.");
}

//////////////////////////////////////////////////////////////////////////////

PrepareDaedalusFiles_Valve::PrepareDaedalusFiles_Valve(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
   addConfigOptionsTo(this);

   _script = "";
   setParameter("Script",&_script);

   _scriptFile = "";
   setParameter("ScriptStartFile",&_scriptFile);

}

//////////////////////////////////////////////////////////////////////////////

void PrepareDaedalusFiles_Valve::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  SimpleMeshAdapterCom::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void PrepareDaedalusFiles_Valve::setup()
{

  SimpleMeshAdapterCom::setup();

}

//////////////////////////////////////////////////////////////////////////////

void PrepareDaedalusFiles_Valve::execute()
{
  CFAUTOTRACE;

  //Write the new daedalus files
  if (PE::GetPE().IsParallel()) {

    PE::GetPE().setBarrier();

    if (PE::GetPE().GetRank () == 0) {
      readWriteFile();
    }

    PE::GetPE().setBarrier();
  }
  else{
    readWriteFile();
  }

}

//////////////////////////////////////////////////////////////////////////////

void PrepareDaedalusFiles_Valve::readWriteFile()
{
  CFAUTOTRACE;

  //Compute current valve lift
  const CFreal currentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
const CFreal _minOpening = 0.001;
const CFreal _maxOpening = 0.012;
const CFreal _rpm = 20000.;

  const CFreal cycleTime = 60./ _rpm;
  const CFreal quarterCycleTime = 0.25 * cycleTime;

  const CFreal currentLift = _minOpening + 0.5 * (_maxOpening - _minOpening) + 0.5 * (_maxOpening - _minOpening) * sin(2.*MathTools::MathConsts::CFrealPi()*(currentTime) / quarterCycleTime);



  const CFreal lift = currentLift;

  boost::filesystem::path newFile =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(_scriptFile);

  SelfRegistPtr<Environment::FileHandlerOutput> newJournalHandle =
     Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = newJournalHandle->open(newFile);

  fout << "Valve_profile_parameters\n";
  fout << "\n";
  fout << "lift " << lift << "\n";
  fout << "\n";
  fout << "VALVE \n";
  fout << "seat_angle      -45.00 \n";
  fout << "seat_lenght       2.1929319102 \n";
  fout << "seat_tangent     -2.00 \n";
  fout << "given_point_X    95.00 \n";
  fout << "given_point_Y     6.00 \n";
  fout << "valve_tangent     0.00 \n";
  fout << "horizontal_part_beginning   115.00 \n";
  fout << "\n";
  fout << "CYLINDER \n";
  fout << "seat_tangent     -0.60 \n";
  fout << "wall_angle       95.00 \n";
  fout << "\n";
  fout << "writing_pts       1 \n";

  newJournalHandle->close();

  Common::OSystem::getInstance().executeCommand(_script);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
