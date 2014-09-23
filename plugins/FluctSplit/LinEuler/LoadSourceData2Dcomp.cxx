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

#include "LoadSourceData2Dcomp.hh"
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

MethodCommandProvider<LoadSourceData2Dcomp,
                     DataProcessingData,
                     LinearizedEulerModule>
aLoadSourceData2DcompProvider("LoadSourceData2Dcomp");

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dcomp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("InputFile","Name of Input File.");
  options.addConfigOption< CFuint >("StartTime","Start of the source injection.");  
  options.addConfigOption< CFreal >("TimeStep","Time-step.");
  options.addConfigOption< CFuint >("nbIter","Number of iterations.");
}

//////////////////////////////////////////////////////////////////////////////

LoadSourceData2Dcomp::LoadSourceData2Dcomp( const std::string& name) :
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

LoadSourceData2Dcomp::~LoadSourceData2Dcomp()
{
}

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dcomp::setup()
{
  DataProcessingCom::setup(); // first call setup of parent class
  
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
LoadSourceData2Dcomp::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_sourceData);

  return result;
}

//////////////////////////////////////////////////////////////////////////////


void LoadSourceData2Dcomp::execute(){
  
  DataHandle<RealVector> sourceData = socket_sourceData.getDataHandle();
  
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
  
  Common::SafePtr< vector<CFuint> > statesIdx = innerTrs->getStatesInTrs();
  
  CFuint nbstates = statesIdx->size();
  
  sourceData.resize(nbstates);
  
  for (CFuint inner=0; inner < nbstates; inner++)
    sourceData[inner].resize(nbEqs+2);
  
//   cout << "Source data size:  " << sourceData[innerNbGeos-1] << "\n";
  
// set the source time
  CFint sourceTime;
  CFreal simulationTime = 100000.*time;
  
  sourceTime = (int)(simulationTime/m_TimeStep);
  
  CFuint nsamples;

  nsamples = sourceTime*m_TimeStep+m_StartTime;

  for(CFuint iCell=0; iCell<nbstates; iCell++) {
       sourceData[iCell][0] = 0.0;
       sourceData[iCell][1] = 0.0;
       sourceData[iCell][2] = 0.0;
       sourceData[iCell][3] = 0.0;
       sourceData[iCell][4] = 0.0;
       sourceData[iCell][5] = 0.0;
    }
  
  std::ostringstream filename;
  filename << m_InputFile << nsamples << ".dat";
  boost::filesystem::path fnameS = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandleS = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  
  ifstream& fileS = fhandleS->open(fnameS);
  
  for(CFuint st=0; st<nbEqs+2; st++) {
    for(CFuint i=0; i<nbstates; i++) {
      
      fileS >> sourceData[i][st];
      }
    }

  fhandleS->close();  

  cout << "Source data loaded for time......  " << nsamples << "\n";
  
  cout << "Source file......  " << fnameS << "\n";

}//execute

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dcomp::unsetup()
{

  DataProcessingCom::unsetup(); // at last call setup of parent class

}

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData2Dcomp::configure (Config::ConfigArgs& args)
{
  DataProcessingCom::configure( args );

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler
    
  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////



