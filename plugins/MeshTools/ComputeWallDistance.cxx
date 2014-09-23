// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "MeshTools/MeshTools.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "MeshTools/ComputeWallDistance.hh"
#include "Common/BadValueException.hh"
#include "Common/NoSuchValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeWallDistance, DataProcessingData, MeshToolsModule> computeWallDistanceProvider("ComputeWallDistance");

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistance::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("OutputFile","Name of Output File to write the distance.");
   options.addConfigOption< std::vector<std::string> >("BoundaryTRS","TRSs to compute the distance to");
}

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistance::ComputeWallDistance(const std::string& name) :
  DataProcessingCom(name),
  socket_wallDistance("wallDistance"),
  // socket_normals("normals"), 
  socket_nodes("nodes"),
  socket_states("states"),
  _tmpVector()
{
  addConfigOptionsTo(this);

  _nameOutputFile = "WallDistance.dat";
  setParameter("OutputFile",&_nameOutputFile);

  _boundaryTRS = vector<std::string>();
  setParameter("BoundaryTRS",&_boundaryTRS);
}

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistance::~ComputeWallDistance()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeWallDistance::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeWallDistance::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  //  result.push_back(&socket_normals); 
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistance::setup()
{
  CFAUTOTRACE;

  _tmpVector.resize(PhysicalModelStack::getActive()->getDim());

  // Get number of states to resize the datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();

  DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
  wallDistance.resize(nbStates);
  wallDistance = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistance::execute()
{
  CFAUTOTRACE;

  CFLog(INFO, "Computing wall distances using simple state-node distance...\n");
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

  const CFuint nbStates = states.size();
  CFreal minimumDistance;
  CFreal distance;
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    minimumDistance = 10.E99;
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size() ; ++iTRS)
    {
      Common::SafePtr<std::vector<CFuint> > nodesInTrs = MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS])->getNodesInTrs();
      for(CFuint iBoundaryNode = 0; iBoundaryNode < nodesInTrs->size() ; ++iBoundaryNode)
      {
//       std::cout << "Comparing node: " << *(states[iState]->getCoordinates()) << " with node: " << *(nodes[(*nodesInTrs)[iBoundaryNode]]) << std::endl;
      _tmpVector = *(nodes[(*nodesInTrs)[iBoundaryNode]]);
      _tmpVector -= states[iState]->getCoordinates();
      distance = _tmpVector.norm2();
      minimumDistance = min(distance, minimumDistance);
      }
    }
    wallDistance[iState] = minimumDistance;
  }

  printToFile();
  CFLog(INFO, "Wall distances computation finished...\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistance::printToFile()
{
  CFAUTOTRACE;

  using namespace boost::filesystem;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

  path file = Environment::DirPaths::getInstance().getResultsDir() / path(_nameOutputFile);
//   path file = Environment::DirPaths::getInstance().getWorkingDir() / path(_nameOutputFile);
  file = Framework::PathAppender::getInstance().appendParallel( file );
  change_extension(file,".dat");

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(file);

  const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();

  fout << "!WALLDISTANCE" << endl;
  fout << "!NBSTATES " << wallDistance.size() << endl;

  if(dim == 2) fout << "!VARIABLES x0 x1 StateID Distance" << endl;
  if(dim == 3) fout << "!VARIABLES x0 x1 x2 StateID Distance" << endl;

  for(CFuint iState=0;iState < wallDistance.size() ; ++iState)
  {
    fout << states[iState]->getCoordinates() << "  " << iState << " " << wallDistance[iState] << endl;
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistance::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  if (_boundaryTRS.size()<1)
    throw Common::BadValueException(FromHere(),"ComputeWallDistance: size of option 'BoundaryTrs' is less then one.");

  Common::SafePtr<MeshData> mesh_data = MeshDataStack::getInstance().getEntryByNamespace(getMethodData().getNamespacePtr());
  std::vector<bool> boundaryTRSfound(_boundaryTRS.size(),false);
  std::vector< std::string > TRSlist=mesh_data->getTRSNameList();

  for(CFuint iTRS = 0; iTRS < _boundaryTRS.size() ; ++iTRS) {
    for(CFuint jTRS = 0; jTRS < TRSlist.size() ; ++jTRS)
      if(_boundaryTRS[iTRS]==TRSlist[jTRS])
        boundaryTRSfound[iTRS]= true;
    if (boundaryTRSfound[iTRS]==false)
      throw Common::NoSuchValueException(FromHere(),"ComputeWallDistance: '"+_boundaryTRS[iTRS]+"'specified for 'BoundaryTrs' option was not found. Maybe typo in Trs name?");
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
