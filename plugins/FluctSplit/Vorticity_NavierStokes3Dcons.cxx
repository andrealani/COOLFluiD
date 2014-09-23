#include <iostream>
#include <string>
#include <fstream>

#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FluctSplit/Vorticity_NavierStokes3Dcons.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

using namespace COOLFluiD::Physics::NavierStokes;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Vorticity_NavierStokes3Dcons,
                      DataProcessingData,
		      FluctSplitSpaceTimeModule>
vorticityns3dconsProvider("Vorticity_NavierStokes3Dcons");

//////////////////////////////////////////////////////////////////////////////

void Vorticity_NavierStokes3Dcons::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("OutputFile","Name of Output File.");
  options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file with aerodynamic coefficients.");
  options.addConfigOption< bool >("AppendTime","Append time to file name.");
  options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");

}

//////////////////////////////////////////////////////////////////////////////

Vorticity_NavierStokes3Dcons::Vorticity_NavierStokes3Dcons( const std::string& name) :
  DataProcessingCom(name),
  m_fhandle(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_vorticity("vorticity"),
  socket_normals("normals"),
  socket_volumes("volumes")
{
  m_fileVorticity = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);
  m_nameOutputFileVorticity = "Vorticity_save.dat";
  setParameter("OutputFile",&m_nameOutputFileVorticity);
 
  m_saveRateVorticity = 1;
  setParameter("SaveRate",&m_saveRateVorticity);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

}

//////////////////////////////////////////////////////////////////////////////

Vorticity_NavierStokes3Dcons::~Vorticity_NavierStokes3Dcons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Vorticity_NavierStokes3Dcons::needsSockets()
{ 
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states); 
  result.push_back(&socket_vorticity);
  result.push_back(&socket_normals);
  result.push_back(&socket_volumes);  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity_NavierStokes3Dcons::setup()
{
  CFAUTOTRACE;
 DataProcessingCom::setup(); // first call setup of parent class

  m_varSet->setup();

}
//////////////////////////////////////////////////////////////////////////////

void Vorticity_NavierStokes3Dcons::execute()
{
  Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();

  const CFuint iter = ssys_status->getNbIter();

    if(!(iter % m_saveRateVorticity))  { computeVorticity(true); }

}


//////////////////////////////////////////////////////////////////////////////
void Vorticity_NavierStokes3Dcons::computeVorticity(bool save){


   Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  const CFreal time = subSysStatus->getCurrentTime();
 DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
 DataHandle<CFreal> vorticity = socket_vorticity.getDataHandle();
 
 DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
 DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();
  
 const CFuint nbstates = states.size();

//     prepare looping over the cells
//     find the inner trs
  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr<TopologicalRegionSet> innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
      break;
    }
  }
  if (innerTrs==0) Common::NoSuchValueException(FromHere(),"Trs with tag 'inner' not found.");
  
    // prepares to loop over cells by getting the GeometricEntityPool
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint innerNbGeos=innerTrs->getLocalNbGeoEnts();
  
  std::vector<CFreal> allVolumes(nbstates);
   
  for (CFuint iState = 0; iState < nbstates; ++iState) {
     vorticity[iState] = 0.0;
     allVolumes[iState] = 0.0;
  }
     
     // compute the gradients over the cells
   for (CFuint iCell=0; iCell<innerNbGeos; iCell++) {
    
//       build the GeometricEntity
      const CFuint cellID=innerTrs->getLocalGeoID(iCell);
      geoData.idx = cellID;
      GeometricEntity& cell = *geoBuilder->buildGE();
      vector<State*> *const statesInCell = cell.getStates();
      const CFuint nbStatesInCell = statesInCell->size();
     
      
      CFreal dUdX = 0.0;
      CFreal dUdY = 0.0;
      CFreal dUdZ = 0.0;
	
      CFreal dVdX = 0.0;
      CFreal dVdY = 0.0;
      CFreal dVdZ = 0.0;
	
      CFreal dWdX = 0.0;
      CFreal dWdY = 0.0;
      CFreal dWdZ = 0.0;
      
   for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
  
        State *const currState = (*statesInCell)[elemState];
	CFuint elemStateID = currState->getLocalID();  
  
        CFreal rho  = (*states[elemStateID])[0];
	CFreal rhou = (*states[elemStateID])[1];
	CFreal rhov = (*states[elemStateID])[2];
	CFreal rhow = (*states[elemStateID])[3];
	CFreal rhoE = (*states[elemStateID])[4];
   
	CFreal u = rhou/rho;
	CFreal v = rhov/rho;
	CFreal w = rhow/rho;

        dUdX +=m_normals[cellID]->getNodalNormComp(elemState,0) *u;
	dUdY +=m_normals[cellID]->getNodalNormComp(elemState,1) *u;
	dUdZ +=m_normals[cellID]->getNodalNormComp(elemState,2) *u;
	
      	dVdX +=m_normals[cellID]->getNodalNormComp(elemState,0) *v;
	dVdY +=m_normals[cellID]->getNodalNormComp(elemState,1) *v;
	dVdZ +=m_normals[cellID]->getNodalNormComp(elemState,2) *v;
	
      	dWdX +=m_normals[cellID]->getNodalNormComp(elemState,0) *w;
	dWdY +=m_normals[cellID]->getNodalNormComp(elemState,1) *w;
	dWdZ +=m_normals[cellID]->getNodalNormComp(elemState,2) *w;
	
   }
   
   CFreal oneover2vol = 1./(2.0*volumes[cellID]);
   dUdX *= oneover2vol;
   dUdY *= oneover2vol;
   dUdZ *= oneover2vol;   
   
   dVdX *= oneover2vol;
   dVdY *= oneover2vol;
   dVdZ *= oneover2vol;      
      
   dWdX *= oneover2vol;
   dWdY *= oneover2vol;
   dWdZ *= oneover2vol;
   
   for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
     State *const currState = (*statesInCell)[elemState];
     CFuint stateID = currState->getLocalID();  
     vorticity[stateID][0] += (dWdY-dVdZ)*volumes[cellID];
     vorticity[stateID][1] += (dUdZ-dWdX)*volumes[cellID];
     vorticity[stateID][2] += (dVdX-dUdY)*volumes[cellID];
     vorticity[stateID][3] += sqrt(vorticity[stateID][0]*vorticity[stateID][0]+vorticity[stateID][1]*vorticity[stateID][1]+vorticity[stateID][2]*vorticity[stateID][2]);
     allVolumes[stateID] += volumes[cellID];
   }

  //release the GeometricEntity
  geoBuilder->releaseGE();
    
  }

  for (CFuint iState = 0; iState < nbstates; ++iState) {
     State *const currState = states[iState];
     CFuint stateIDGlobal = currState->getLocalID();
     vorticity[stateIDGlobal] /= allVolumes[stateIDGlobal];
  }

}
  

//////////////////////////////////////////////////////////////////////////////

void Vorticity_NavierStokes3Dcons::unsetup()
{ 

}

//////////////////////////////////////////////////////////////////////////////
void Vorticity_NavierStokes3Dcons::prepareOutputFileVorticity()
{

  using boost::filesystem::path;

  cf_assert (!m_fileVorticity->isopen());
  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileVorticity);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime);


  ofstream& fout = m_fileVorticity->open(file);

  fout << "TITLE  =  Vorticity in the whole field" << "\n";
  fout << "VARIABLES = x y z qcrit" << "\n";
  
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity_NavierStokes3Dcons::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);

  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  std::string varSetName = "Euler3DCons";
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(m_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

    } //namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

