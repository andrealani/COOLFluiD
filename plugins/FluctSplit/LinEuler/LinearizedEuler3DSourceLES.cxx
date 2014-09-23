#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <numeric>
#include <algorithm>

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"

#include "LinearizedEuler3DSourceLES.hh"
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

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler3DSourceLES,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler3DSTProviderLES("LinEuler3DSourceLES");
//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceLES::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("InputFile","Name of Input File.");
  options.addConfigOption< CFuint >("StartTime","Start of the source injection.");  
  options.addConfigOption< CFreal >("TimeStep","Time-step.");
  options.addConfigOption< CFuint >("nbIter","Number of iterations.");
  options.addConfigOption< bool >("ComputeSources","Compute the sources from LES/DNS?");
  options.addConfigOption< bool >("Interpolation","Interpolation of meshes.");
  options.addConfigOption< CFreal >("MeanDensity","Meanflow density (constant for incompressible flows.)");
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler3DSourceLES::LinearizedEuler3DSourceLES(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_states("states"),
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
  
  m_ComputeSources = false;
  setParameter("ComputeSources",&m_ComputeSources);

  m_Interpolation = false;
  setParameter("Interpolation",&m_Interpolation);
  
  m_meanDensity = 0.0;
  setParameter("MeanDensity",&m_meanDensity);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedEuler3DSourceLES::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);
  result.push_back(&socket_normals);
  result.push_back(&socket_sourceData);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler3DSourceLES::~LinearizedEuler3DSourceLES()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceLES::setup()
{
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());

    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();
//   DistributionData& ddata = getMethodData().getDistributionData();

  // prepares to loop over cells by getting the GeometricEntityPool
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbstates=states.size();
  
  const CFreal rho0 = m_meanDensity;
  
  // prepare looping over the cells
    // find the inner trs
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
  
  // the sources defined over an element
  std::vector< RealVector > sourceData;
  sourceData.resize(innerNbGeos);
  for(CFuint i=0; i<innerNbGeos; i++)
  sourceData[i].resize(nbEqs);
  
// Initialize the sourceData with zeros

  CFuint nbstatesLES = 0;
  RealVector CTable(nbstates);

  if(m_Interpolation) {
    
 // read the CTable 

  boost::filesystem::path fCT = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path("CTable.dat");
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

  ifstream& fileCT = fhandle->open(fCT);
  
  CFuint number;
  
  for(CFuint ct=0; ct<nbstates; ct++) {
    fileCT >> number;
    fileCT >> CTable[ct];
  }  

  fhandle->close();  

  cout << "\n The connectivity table is loaded. \n";
  
  }
  
if(m_ComputeSources) {
  
  std::vector< RealVector > meanvel, meanstress;
  std::vector< RealVector > meanvelLEE, meanstressLEE;
  
  meanvelLEE.resize(nbstates);
  meanstressLEE.resize(nbstates);

  for (CFuint i = 0; i < nbstates; ++i) {
    meanvelLEE[i].resize(nbdim);
    meanstressLEE[i].resize(3*nbdim-3);
  }

// read the mean values    

  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path("sourcemeans.dat");
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

  ifstream& file = fhandle->open(fname);

  std::string dummy;
  char ch;
  CFuint NNode = 0;

/* header */
  getline (file,dummy);


  for(int i=1; i< 10; i++) {
  	getline (file,dummy);
  }

  getline (file,dummy);
  getline (file,dummy);
  
  file >> ch >> ch; 
  file >> NNode;
  
  nbstatesLES = NNode;
  
  meanvel.resize(NNode);
  meanstress.resize(NNode);

  for (CFuint i = 0; i < NNode; ++i) {
    meanvel[i].resize(nbdim);
    meanstress[i].resize(3*nbdim-3);
  }

  for(CFint i=1; i<4; i++)
    getline (file,dummy);

/* velocities */
  for(CFuint dim=0; dim<nbdim; dim++) {
    for(CFuint i=0; i<NNode; i++)
      file >> meanvel[i][dim];
  }
  
  for(CFuint st=0; st<3*nbdim-3; st++) {
    for(CFuint i=0; i<NNode; i++)
      file >> meanstress[i][st];
    }

  fhandle->close();
  
  cout << "\n The mean quantities are red for the source reconstruction. \n";
  
  
// interpolate if neccessary
  if(m_Interpolation) {
    for (CFuint i = 0; i < nbstates; ++i) {
      meanvelLEE[i] = meanvel[CTable[i]];
      meanstressLEE[i] = meanstress[CTable[i]];
    }
  }
  else {
    for (CFuint i = 0; i < nbstates; ++i) {
      meanvelLEE[i] = meanvel[i];
      meanstressLEE[i] = meanstress[i];
    }
  }
  

// read the instantaneous velocities


  std::vector< RealVector > instVEL;
  instVEL.resize(nbstatesLES);
  for(CFuint i=0; i<nbstatesLES; i++)
    instVEL[i].resize(nbdim);
  
  std::vector< RealVector > instVELLEE;
  instVELLEE.resize(nbstates);
  for(CFuint i=0; i<nbstates; i++)
    instVELLEE[i].resize(nbdim);
  
// loop over the time-instances

  for(CFuint iIter=0; iIter<m_nbIter; iIter++) {
  
  CFuint nsamples = m_StartTime+iIter*m_TimeStep;
  
  cout << "Reconstructing the sources for time:........... " << nsamples << "\n";
  
  FILE* fileV;
  
  char fnameV[256];
  
  sprintf(fnameV,"velo_%06d.dat",nsamples);
  
  
  if( (fileV = fopen( fnameV, "rb")) == NULL) {
   cout << "\nCannot open velocity file!!! \n";
  }
  else {
  
    double *buffer;
    buffer=(double*)calloc(1,sizeof(double));

    for(CFuint nbs=0; nbs<nbstatesLES; nbs++) {
      for(CFuint dim=0; dim<nbdim; dim++) {
        fread(buffer,sizeof(double),1,fileV);
        instVEL[nbs][dim] = *buffer;
      }
    }
    fclose(fileV);
    free(buffer);
  }
  
// interpolate if neccessary
  if(m_Interpolation) {
    for (CFuint i = 0; i < nbstates; ++i) {
      instVELLEE[i] = instVEL[CTable[i]];
    }
  }
  else {
    for (CFuint i = 0; i < nbstates; ++i) {
      instVELLEE[i] = instVEL[i];
    }
  }

// compute the perturbation velocities    
  
  std::vector< RealVector > perturbVEL;
  perturbVEL.resize(nbstates);
  for(CFuint i=0; i<nbstates; i++)
    perturbVEL[i].resize(nbdim);
  
  for(CFuint nbs=0; nbs<nbstates; nbs++) {
    for(CFuint dim=0; dim<nbdim; dim++) { 
      perturbVEL[nbs][dim] = instVELLEE[nbs][dim]-meanvelLEE[nbs][dim];
    }
  }

// compute the gradients and the sources over the cells
  for (CFuint iCell=0; iCell<innerNbGeos; iCell++) {
    
  // build the GeometricEntity
    const CFuint cellID=innerTrs->getLocalGeoID(iCell);
    geoData.idx = cellID;
    GeometricEntity& cell = *geoBuilder->buildGE();
    vector<State*> *const statesInCell = cell.getStates();
    const CFuint nbStatesInCell = statesInCell->size();

 // compute the meanstress derivatives    
   CFreal drhouudx_mean = 0.0;
   CFreal drhovvdy_mean = 0.0;
   CFreal drhowwdz_mean = 0.0;
   
   CFreal drhouvdx_mean = 0.0;
   CFreal drhouvdy_mean = 0.0;
   
   CFreal drhouwdx_mean = 0.0;
   CFreal drhouwdz_mean = 0.0;
   
   CFreal drhovwdy_mean = 0.0;
   CFreal drhovwdz_mean = 0.0;
      
   for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
      State *const currState = (*statesInCell)[elemState];
      CFuint stateID = currState->getLocalID();
     
      	drhouudx_mean +=m_normals[cellID]->getNodalNormComp(elemState,0) *meanstressLEE[stateID][0];
	drhovvdy_mean +=m_normals[cellID]->getNodalNormComp(elemState,1)*meanstressLEE[stateID][3];
	drhowwdz_mean +=m_normals[cellID]->getNodalNormComp(elemState,2)*meanstressLEE[stateID][5];
	
	drhouvdx_mean +=m_normals[cellID]->getNodalNormComp(elemState,0)*meanstressLEE[stateID][1];
	drhouvdy_mean +=m_normals[cellID]->getNodalNormComp(elemState,1)*meanstressLEE[stateID][1];
	
	drhouwdx_mean +=m_normals[cellID]->getNodalNormComp(elemState,0)*meanstressLEE[stateID][2];
	drhouwdz_mean +=m_normals[cellID]->getNodalNormComp(elemState,2)*meanstressLEE[stateID][2];
	
	drhovwdy_mean +=m_normals[cellID]->getNodalNormComp(elemState,1)*meanstressLEE[stateID][4];
	drhovwdz_mean +=m_normals[cellID]->getNodalNormComp(elemState,2)*meanstressLEE[stateID][4];
   }
   
   CFreal oneover2vol = 1./(2.0*volumes[cellID]);
   drhouudx_mean *= oneover2vol;
   drhovvdy_mean *= oneover2vol;
   drhowwdz_mean *= oneover2vol;
   
   drhouvdx_mean *= oneover2vol;
   drhouvdy_mean *= oneover2vol;

   drhouwdx_mean *= oneover2vol;
   drhouwdz_mean *= oneover2vol;
   
   drhovwdy_mean *= oneover2vol;
   drhovwdz_mean *= oneover2vol;
   
   
//compute the stresses and their derivatives
   CFreal rhouu = 0.0;
   CFreal rhovv = 0.0;
   CFreal rhoww = 0.0;
   
   CFreal rhouv = 0.0;
   CFreal rhouw = 0.0;
   CFreal rhovw = 0.0;
   
   CFreal drhouudx = 0.0;
   CFreal drhovvdy = 0.0;
   CFreal drhowwdz = 0.0;
   
   CFreal drhouvdx = 0.0;
   CFreal drhouvdy = 0.0;
   
   CFreal drhouwdx = 0.0;
   CFreal drhouwdz = 0.0;
   
   CFreal drhovwdy = 0.0;
   CFreal drhovwdz = 0.0;
   
   for (CFuint iState = 0; iState <nbStatesInCell; ++iState) {
        State *const currState = (*statesInCell)[iState];
	CFuint stateID = currState->getLocalID();
  
	rhouu = rho0*perturbVEL[stateID][0]*perturbVEL[stateID][0];
	rhovv = rho0*perturbVEL[stateID][1]*perturbVEL[stateID][1];
	rhoww = rho0*perturbVEL[stateID][2]*perturbVEL[stateID][2];
	
	rhouv = rho0*perturbVEL[stateID][0]*perturbVEL[stateID][1];
	rhouw = rho0*perturbVEL[stateID][0]*perturbVEL[stateID][2];
	rhovw = rho0*perturbVEL[stateID][1]*perturbVEL[stateID][2];
     
       	drhouudx +=m_normals[cellID]->getNodalNormComp(iState,0)*rhouu;
	drhovvdy +=m_normals[cellID]->getNodalNormComp(iState,1)*rhovv;
	drhowwdz +=m_normals[cellID]->getNodalNormComp(iState,2)*rhoww;
	
	drhouvdx +=m_normals[cellID]->getNodalNormComp(iState,0)*rhouv;
	drhouvdy +=m_normals[cellID]->getNodalNormComp(iState,1)*rhouv;
	
	drhouwdx +=m_normals[cellID]->getNodalNormComp(iState,0)*rhouw;
	drhouwdz +=m_normals[cellID]->getNodalNormComp(iState,2)*rhouw;
	
	drhovwdy +=m_normals[cellID]->getNodalNormComp(iState,1)*rhovw;
	drhovwdz +=m_normals[cellID]->getNodalNormComp(iState,2)*rhovw;
   }
   
   drhouudx *= oneover2vol;
   drhovvdy *= oneover2vol;
   drhowwdz *= oneover2vol;
   
   drhouvdx *= oneover2vol;
   drhouvdy *= oneover2vol;

   drhouwdx *= oneover2vol;
   drhouwdz *= oneover2vol;
   
   drhovwdy *= oneover2vol;
   drhovwdz *= oneover2vol;
   
   
// compute the source terms
  sourceData[cellID][0] = 0.0;
  sourceData[cellID][1] = -(drhouudx+drhouvdy+drhouwdz) + drhouudx_mean+drhouvdy_mean+drhouwdz_mean;
  sourceData[cellID][2] = -(drhouvdx+drhovvdy+drhovwdz) + drhouvdx_mean+drhovvdy_mean+drhovwdz_mean;
  sourceData[cellID][3] = -(drhouwdx+drhovwdy+drhowwdz) + drhouwdx_mean+drhovwdy_mean+drhowwdz_mean;
  sourceData[cellID][4] = 0.0;

    geoBuilder->releaseGE();

  }
  
// save source data to file
  
  std::ostringstream filename;
  filename << m_InputFile << nsamples << ".dat";
  boost::filesystem::path fnameS = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleS = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  
  ofstream& fileS = fhandleS->open(fnameS);
  

//write just the sources in block structure
  for(CFuint st=1; st<nbEqs; st++) {
    for(CFuint i=0; i<innerNbGeos; i++)
      fileS << sourceData[i][st] << "\t";
    fileS << "\n";
    }

  fhandleS->close();  
    
  }
  
}

cout << "\n The reconstructed sources are saved. \n";  
  
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceLES::setVarSet(SafePtr<ConvectiveVarSet> varSet)
{
  _varSet = varSet.d_castTo<LinEuler2DVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceLES::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData&m_normals)
{
  
  DataHandle<RealVector> sourceData = socket_sourceData.getDataHandle();
  
  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint cellID = ddata.cellID;
  
  source = sourceData[cellID];
  
  
  

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
