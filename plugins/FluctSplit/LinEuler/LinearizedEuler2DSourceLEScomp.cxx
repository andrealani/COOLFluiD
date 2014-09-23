
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <numeric>
#include <algorithm>

#include <string>
#include <fstream>
#include <iostream>

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"

#include "LinearizedEuler2DSourceLEScomp.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluctSplitLinEuler.hh"
#include "FluctSplit/FluctuationSplitData.hh"

#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"
#include "Framework/PathAppender.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler2DSourceLEScomp,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler2DSTProviderLEScomp("LinEuler2DSourceLEScomp");
//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceLEScomp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("InputFile","Name of Input File.");
  options.addConfigOption< std::string >("prePath","Name of the path to sources.");
  options.addConfigOption< std::string >("postPath","Name of the sources files.");
  options.addConfigOption< std::string >("postMeanPath","Name of the mean sources files.");
  options.addConfigOption< CFreal >("StartTime","Start of the source injection.");
  options.addConfigOption< CFreal >("TimeStep","Time-step.");
  options.addConfigOption< CFreal >("radius","Source effectiveness radius");
  options.addConfigOption< CFuint >("nbIter","Number of iterations.");
  options.addConfigOption< bool >("ComputeSources","Compute the sources from LEScomp/DNS?");
  options.addConfigOption< bool >("Interpolation","Interpolation of meshes.");
  options.addConfigOption< bool >("CTable","whether to compute the CTable or read it");
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceLEScomp::LinearizedEuler2DSourceLEScomp(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_states("states"),
  socket_sourceSave("sourceSave")
{
  addConfigOptionsTo(this);

  m_InputFile = "sourceData";
  setParameter("InputFile",&m_InputFile);

  m_prePath = "./";
  setParameter("prePath",&m_prePath);

  m_postPath = "/Si_constantPlane.raw";
  setParameter("postPath",&m_postPath);

  m_postMeanPath = "/SiMean_constantPlane.raw";
  setParameter("postMeanPath",&m_postMeanPath);

  m_StartTime = 1;
  setParameter("StartTime",&m_StartTime);
  
  m_radius = 0.05;
  setParameter("radius",&m_radius);
  
  m_TimeStep = 1;
  setParameter("TimeStep",&m_TimeStep);
  
  m_nbIter = 1;
  setParameter("nbIter",&m_nbIter);
  
  m_ComputeSources = false;
  setParameter("ComputeSources",&m_ComputeSources);

  m_Interpolation = false;
  setParameter("Interpolation",&m_Interpolation);

  m_CTable = false;
  setParameter("CTable",&m_CTable);
  
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedEuler2DSourceLEScomp::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);
  result.push_back(&socket_normals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////
std::vector<Common::SafePtr<BaseDataSocketSource> >
LinearizedEuler2DSourceLEScomp::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_sourceSave);
  return result;
}
//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceLEScomp::~LinearizedEuler2DSourceLEScomp()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceLEScomp::setup()
{

  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();

  DataHandle<RealVector> sourceSave = socket_sourceSave.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();
  
  //past time
  m_pastTime = -1.0;

  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbstates=states.size();

  sourceData.resize(nbstates);
  for(CFuint i=0; i<nbstates; i++) {
     sourceData[i].resize(nbdim*2);
  }
  
// Initialize the sourceData with zeros
  for(CFuint i=0; i<nbstates; i++)
      for(CFuint j=0; j<nbdim*2; j++)
	sourceData[i][j]=0.0;

  const CFuint nbcells = MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  //flag for the mean reading part
  readMean_flag = true;
  CTable.resize(nbstates);
  nbstatesLEScomp = 0;
  
  //read mean stress
  double DUMMYdbl;//to avoid storing x,y,z from extracted data
  string dummy;
  ostringstream ost;

 //Computes the last time in order to get the mean sources
 ost.str("");
 ost.clear();

 CFreal lastTime = m_StartTime+(m_nbIter-1)*m_TimeStep;

 ost << lastTime;
 string fullPath = m_prePath + ost.str() + m_postMeanPath;  
 
 boost::filesystem::path fSi = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(fullPath);
 
 if(readMean_flag){
  readMean_flag = false;
 
 }

 Common::SelfRegistPtr<Environment::FileHandlerInput> fhandleSmean = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

 ifstream& fileSiMean = fhandleSmean->open(fSi);

// header 
  fileSiMean >> dummy >> dummy >> dummy >> nbstatesLEScomp;
  getline(fileSiMean,dummy);
  getline(fileSiMean,dummy);
  
  std::vector< RealVector > LEScoord, LESsource;
  LEScoord.resize(nbstatesLEScomp);
  LESsource.resize(nbstatesLEScomp);  
  for (CFuint i = 0; i < nbstatesLEScomp; ++i) {
     LEScoord[i].resize(3);
     LESsource[i].resize(3);     
  }

 for(CFuint nbs=0; nbs<nbstatesLEScomp; nbs++) {
   LEScoord[nbs][0] = nbs;
   fileSiMean >> LEScoord[nbs][1] >> LEScoord[nbs][2] >> DUMMYdbl >> LESsource[nbs][0] >> LESsource[nbs][1];
   getline(fileSiMean,dummy);
 }

 fhandleSmean->close();  
     
 if(m_Interpolation) {
  if(m_CTable){
    
    ifstream inFileCT;
    
    inFileCT.open("CTable.dat");
    cout << "Reading existing Connectivity table \n";

    for(CFuint ct=0; ct<nbstates; ct++) {
     inFileCT >> CTable[ct];
     getline(inFileCT,dummy);
    }
 
    inFileCT.close();
  }
  else
  {
    std::vector< RealVector > LEEcoord;
    LEEcoord.resize(nbstates);
    for (CFuint i = 0; i < nbstates; ++i) {
      LEEcoord[i].resize(3);
      }
  //get the LEE mesh coordinate
    for (CFuint iState = 0; iState < nbstates; ++iState)  {
      const CFreal x = (states[iState]->getCoordinates())[XX];
      const CFreal y = (states[iState]->getCoordinates())[YY];
   
      LEEcoord[iState][0]=iState;
      LEEcoord[iState][1]=x;
      LEEcoord[iState][2]=y;
    }
 
 //************* compute the Connectivity table and write it ******************************
    
    boost::filesystem::path fCT = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path("CTable.dat");

    Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleOutCT = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fileCT = fhandleOutCT->open(fCT);

    for(CFuint ctlee=0; ctlee<nbstates; ctlee++) {
      CFreal dist = 1.0e+06; 
      for(CFuint ctles=0; ctles<nbstatesLEScomp; ctles++) {
        if(sqrt( pow((LEEcoord[ctlee][1]-LEScoord[ctles][1]),2) + pow((LEEcoord[ctlee][2]-LEScoord[ctles][2]),2) ) < dist){
          dist = sqrt( pow((LEEcoord[ctlee][1]-LEScoord[ctles][1]),2) + pow((LEEcoord[ctlee][2]-LEScoord[ctles][2]),2) );
          CTable[ctlee]=LEScoord[ctles][0];
        }
      }
      fileCT << CTable[ctlee] << endl;
    }
    fhandleOutCT->close();  
   }
   cout << "\n The connectivity table is ready. \n";
   for (CFuint i = 0; i < nbstates; ++i) {
	for(CFuint dim=0; dim<2; dim++) {
		sourceData[i][dim] = LESsource[(CTable[i])][dim] ;
	}
   }  
  }
 else {
  for (CFuint i = 0; i < nbstates; ++i) {
	for(CFuint dim=0; dim<2; dim++) {
		sourceData[i][dim] = LESsource[i][dim] ;
	}   
  }
 }

}//void Linearized...comp::setup

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceLEScomp::readOFSource()  //function which reads the sources and stores them in sourceData
{
    
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbstates=states.size();
  
  std::vector< RealVector > Si;
  Si.resize(nbstatesLEScomp);
  for(CFuint i=0; i<nbstatesLEScomp; i++)
    Si[i].resize(nbdim*2);

  for(CFuint i=0; i<nbstatesLEScomp; i++){
      for(CFuint j=0; j<nbdim*2; j++){
	Si[i][j]=0.0;
      }
  }

 DistributionData& ddata = getMethodData().getDistributionData();
 const CFreal currentTime = ddata.time;
 string fnameS;
    
 //clear the ostringstream
 ostringstream ost;
 ost.str("");
 ost.clear();

 //clear the string
 fnameS.clear();
 ost << currentTime + m_StartTime;
 fnameS = ost.str(); 
 
 //construct the full path:
 string fullPath = m_prePath + fnameS + m_postPath;

  //loading stress
 
 boost::filesystem::path fSi = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(fullPath); 

  //boost::filesystem::path 
  fSi = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(fullPath);
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandleS = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  
  string dummy;
  CFreal DUMMYdbl;

  ifstream&  fileSi = fhandleS->open(fSi);

 for(int i=0; i<2; i++) {
  getline(fileSi,dummy);
 }
 for(CFuint nbs=0; nbs<nbstatesLEScomp; nbs++) {
   fileSi >> DUMMYdbl >> DUMMYdbl >> DUMMYdbl >> Si[nbs][2] >> Si[nbs][3];
   getline(fileSi,dummy);
 }

 fhandleS->close();

// interpolate if neccessary
 if(m_Interpolation) {

   for (CFuint i = 0; i < nbstates; ++i) {
	for(CFuint dim=2; dim<nbdim*2; dim++) {
		sourceData[i][dim] = Si[(CTable[i])][dim] ;
		}
	}
  }
  else {
	for (CFuint i = 0; i < nbstates; ++i) {
		for(CFuint dim=2; dim<nbdim*2; dim++) {
			sourceData[i][dim] = Si[i][dim] ;
		}
	}
  }

// // save source data to file
//   
//   std::ostringstream filename;
//   filename << "LoadedSources" << currentTime << ".dat";
//   boost::filesystem::path fnameSL = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
//   
//   Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleSL = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
//   
//   ofstream& fileSL= fhandleSL->open(fnameSL);
//   
// 
// //write just the sources in block structure
//   for(CFuint st=0; st<nbstates; st++) {
//     CFreal x = (states[st]->getCoordinates())[XX];
//     CFreal y = (states[st]->getCoordinates())[YY];
//     fileSL << x << "\t" << y << "\t";
//     
//     for(CFuint i=0; i<4; i++) {
//       fileSL << sourceData[st][i] << "\t";
//     }
//     fileSL << "\n";
//    }
// 
//   fhandleSL->close();
   
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceLEScomp::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData&m_normals)
{
    
  DistributionData& ddata = getMethodData().getDistributionData();
  //current time of the computation 
  
  const CFreal currentTime = ddata.time;
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;
  
  //to check if we're at the same time step ; if yes, no need to read the source then. Need to add a condition for source data last time?
  if (currentTime > m_pastTime){
  readOFSource();
  m_pastTime = currentTime;
  }
  const CFuint cellID = ddata.cellID;
   
  const CFuint nbStatesInCell = ddata.states->size();
  vector<State*> *const statesInCell = cell->getStates();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle(); 
  vector<State*>& states = *ddata.states;
  
  
  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;
  
  CFreal radius = 0.0;
      
   
   for (CFuint iState = 0; iState <nbStatesInCell; ++iState) {
        State *const currState = (*statesInCell)[iState];
	CFuint stateID = currState->getLocalID();
        const CFreal x = (states[iState]->getCoordinates())[XX];
        const CFreal y = (states[iState]->getCoordinates())[YY];

	radius = sqrt(x*x+y*y);

	//for now, the average is taken, a better interpolation should be done here
	if (radius < m_radius) {
          source[1] += sourceData[stateID][2]-sourceData[stateID][0];
          source[2] += sourceData[stateID][3]-sourceData[stateID][1];
	}
   }
  
   const CFuint dimSource = source.size();
  
   for (CFuint iSource = 0; iSource < dimSource; ++iSource) {
     source[iSource] /= nbStatesInCell;
     source[iSource] *= volumes[cellID]*dt;
   }

}//void computesourceFSM

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
