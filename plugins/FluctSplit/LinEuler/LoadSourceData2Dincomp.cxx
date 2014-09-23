#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <numeric>
#include <algorithm>

#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/CFL.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerInput.hh"

#include "LoadSourceData2Dincomp.hh"
#include "LinEuler/LinearizedEuler.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
 
    namespace Physics {
      
      namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LoadSourceData2Dincomp,
                     DataProcessingData,
                     LinearizedEulerModule>
aLoadSourceData2DincompProvider("LoadSourceData2Dincomp");

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dincomp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("InputFile","Name of Input File.");
  options.addConfigOption< CFuint >("StartTime","Start of the source injection.");  
  options.addConfigOption< CFreal >("TimeStep","Time-step.");
  options.addConfigOption< CFuint >("nbIter","Number of iterations.");
}

//////////////////////////////////////////////////////////////////////////////

LoadSourceData2Dincomp::LoadSourceData2Dincomp( const std::string& name) :
  DataProcessingCom(name),
  socket_sourceData("sourceData")
{
  addConfigOptionsTo(this);

  m_InputFile = "sourceData";
  setParameter("InputFile",&m_InputFile);

  m_StartTime = 1;
  setParameter("StartTime",&m_StartTime);
  
  m_TimeStep = 1;
  setParameter("TimeStep",&m_TimeStep);
  
  m_nbIter = 1;
  setParameter("nbIter",&m_nbIter);
  
}

//////////////////////////////////////////////////////////////////////////////

LoadSourceData2Dincomp::~LoadSourceData2Dincomp()
{
}

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dincomp::setup()
{
  DataProcessingCom::setup(); // first call setup of parent class
  
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
LoadSourceData2Dincomp::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_sourceData);

  return result;
}

//////////////////////////////////////////////////////////////////////////////


void LoadSourceData2Dincomp::execute(){
  
  DataHandle<RealVector> sourceData = socket_sourceData.getDataHandle();
  
// read source data from file

// find the file needed

  const CFreal time = SubSystemStatusStack::getActive()->getCurrentTime();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr<TopologicalRegionSet> innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
      break;
    }
  }
  if (innerTrs==0) Common::NoSuchValueException(FromHere(),"Trs with tag 'inner' not found.");
  
  const CFuint innerNbGeos=innerTrs->getLocalNbGeoEnts();
  
  sourceData.resize(innerNbGeos);
  
  for (CFuint inner=0; inner < innerNbGeos; inner++)
    sourceData[inner].resize(nbEqs);
  
//   cout << "Source data size:  " << sourceData[innerNbGeos-1] << "\n";
  
// set the source time
  CFint sourceTime;
  CFreal simulationTime = 100000.*time;
  
  sourceTime = (int)(simulationTime/m_TimeStep);
  
  CFuint nsamples;

  nsamples = sourceTime+m_StartTime;

//   cout << "Simulation time:  " << simulationTime  <<  "   Source time:  " << nsamples;
  
//   std::ostringstream filename;
//   filename << m_InputFile << nsamples << ".dat";

  for(CFuint iCell=0; iCell<innerNbGeos; iCell++) {
       sourceData[iCell][0] = 0.0;
    } 
  

  char filename[256];
  sprintf(filename,"source_%d.dat",nsamples);
  
  
//   cout << "File to open: " << filename << "\n";
   
  FILE* fileS;
    
  if( (fileS = fopen( filename, "r")) == NULL) {
   cout << "\nCannot open source file!!! \n";
  }
  else {
   for(CFuint st=1; st<nbEqs; st++) {
    for(CFuint i=0; i<innerNbGeos; i++)
      
      fread(&sourceData[i][st], sizeof(double), 1, fileS);
    }
    
      cout << "Source data loaded from....  " << filename << "\n";
  
    fclose(fileS);
  }

}//execute

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dincomp::unsetup()
{

  DataProcessingCom::unsetup(); // at last call setup of parent class

}

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dincomp::configure (Config::ConfigArgs& args)
{
  DataProcessingCom::configure( args );

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler
    
  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////



