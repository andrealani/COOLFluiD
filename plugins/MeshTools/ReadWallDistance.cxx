// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"

#include "MeshTools/MeshTools.hh"
#include "MeshTools/ReadWallDistance.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ReadWallDistance, DataProcessingData, MeshToolsModule> readWallDistanceProvider("ReadWallDistance");

//////////////////////////////////////////////////////////////////////////////

void ReadWallDistance::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("InputFile","Name of Input File to read the distance.");
}

//////////////////////////////////////////////////////////////////////////////

ReadWallDistance::ReadWallDistance(const std::string& name) :
  DataProcessingCom(name),
  socket_wallDistance("wallDistance"),
  socket_states("states")
{
  addConfigOptionsTo(this);

  _nameInputFile = "WallDistance.dat";
  setParameter("InputFile",&_nameInputFile);
}

//////////////////////////////////////////////////////////////////////////////

ReadWallDistance::~ReadWallDistance()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ReadWallDistance::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ReadWallDistance::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ReadWallDistance::setup()
{
  CFAUTOTRACE;

  // Get number of states to resize the datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();

  DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
  wallDistance.resize(nbStates);
  wallDistance = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void ReadWallDistance::execute()
{
  CFAUTOTRACE;

  using namespace boost::filesystem;

  // Read the file and fill in the wallDistance datahandle
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

  const CFuint nbStates = wallDistance.size();

  path file = Environment::DirPaths::getInstance().getWorkingDir() / path(_nameInputFile);
  file = Framework::PathAppender::getInstance().appendParallel( file );
#ifdef CF_HAVE_BOOST_1_85
  file.replace_extension(".dat");
#else
   change_extension(file,".dat");
#endif

  SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(file);

  std::string line;
  vector<std::string> words;

  //Check content of first line
  getline(fin,line);
  words = Common::StringOps::getWords(line);
  if(words.size() != 1) {
    throw BadFormatException (FromHere(),"Wrong number of parameters in 1st line of file: " + file.string());
  }

  // Check number of nodes identifier
  if(words[0] != "!WALLDISTANCE"){
    throw BadFormatException (FromHere(),"Expecting !WALLDISTANCE identifier in 2nd line of " + file.string());
  }

  //Check content of second line
  getline(fin,line);
  words = Common::StringOps::getWords(line);
  if(words.size() != 2)  {
    throw BadFormatException (FromHere(),"Wrong number of parameters in 2nd line of file: " + file.string());
  }

  // Check number of states identifier
  if(words[0] != "!NBSTATES"){
    throw BadFormatException (FromHere(),"Expecting !NBSTATES identifier in 2nd line of " + file.string());
  }

  // Check agreement of the number of states
  CFuint nbStatesRead = Common::StringOps::from_str<CFint>(words[1]);
  if(nbStatesRead != nbStates){
    throw BadFormatException (FromHere(),"Number of states in file " + file.string() + " differs from number of number of states in mesh");
  }

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  //Check content of third line
  getline(fin,line);
  words = Common::StringOps::getWords(line);

  CFreal temp;
  CFuint ID;

  if( nbDim == 2 ){
    if(words.size() != 5) {
      throw BadFormatException (FromHere(),"Wrong number of parameters in 3rd line of file: " + file.string() + "   Surely in 2D?") ;
    }

    if(words[0] != "!VARIABLES") throw BadFormatException (FromHere(),"Expecting !VARIABLES identifier in 3rd line of " + file.string());
    if(words[1] != "x0")         throw BadFormatException (FromHere(),"Expecting 'x0' in 3rd line of " + file.string());
    if(words[2] != "x1")         throw BadFormatException (FromHere(),"Expecting 'x1' in 3rd line of " + file.string());
    if(words[3] != "StateID")    throw BadFormatException (FromHere(),"Expecting 'StateID' identifier in 3rd line of " + file.string());
    if(words[4] != "Distance")   throw BadFormatException (FromHere(),"Expecting 'Distance' identifier in 3rd line of " + file.string());

    //Read the wall distances
    for (CFuint i = 0; i < nbStates; ++i) {
      // due to new format of distance file some arguments are redundant
      fin >> temp >> temp >> ID; // coordinate x0 -- coordinate x1 -- ID of state
      fin >> wallDistance[ID];
    }
  }


  if( nbDim == 3 ){
    if(words.size() != 6) {
      throw BadFormatException (FromHere(),"Wrong number of parameters in 3rd line of file: " + file.string() + "   Surely in 3D?") ;
    }

    if(words[0] != "!VARIABLES") throw BadFormatException (FromHere(),"Expecting !VARIABLES identifier in 3rd line of " + file.string());
    if(words[1] != "x0")         throw BadFormatException (FromHere(),"Expecting 'x0' in 3rd line of " + file.string());
    if(words[2] != "x1")         throw BadFormatException (FromHere(),"Expecting 'x1' in 3rd line of " + file.string());
    if(words[3] != "x2")         throw BadFormatException (FromHere(),"Expecting 'x2' in 3rd line of " + file.string());
    if(words[4] != "StateID")    throw BadFormatException (FromHere(),"Expecting 'StateID' identifier in 3rd line of " + file.string());
    if(words[5] != "Distance")   throw BadFormatException (FromHere(),"Expecting 'Distance' identifier in 3rd line of " + file.string());

    //Read the wall distances
    for (CFuint i = 0; i < nbStates; ++i) {
      fin >> temp >> temp >> temp >> ID; // coordinate x0 -- coordinate x1 -- coordinate x2 -- ID of state
      fin >> wallDistance[ID];
    }
  }
  
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void ReadWallDistance::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ReadWallDistance::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
