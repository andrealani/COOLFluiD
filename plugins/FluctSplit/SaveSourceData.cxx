#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FluctSplit/SaveSourceData.hh"

#include "NavierStokes/Euler2DCons.hh"
#include "NavierStokes/Euler3DCons.hh"

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

MethodCommandProvider<SaveSourceData,
                      DataProcessingData,
                      FluctSplitModule>
aSaveSourceDataProvider("SaveSourceData");
//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("AppendTime","Append time to file name.");
  options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");  
  options.addConfigOption< CFuint >("SaveRateMean","Compute and save mean stress gradients.");
  options.addConfigOption< bool >("WriteCoordinates","Write coordinates.dat file.");
}

//////////////////////////////////////////////////////////////////////////////

SaveSourceData::SaveSourceData( const std::string& name) :
  DataProcessingCom(name),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_nodes("nodes"),
  socket_states("states")

{
  addConfigOptionsTo(this);
  
  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
  
  m_saveratemean = 1;
  setParameter("SaveRateMean",&m_saveratemean);

  m_writecoordinates = false;
  setParameter("WriteCoordinates",&m_writecoordinates);
  
}

//////////////////////////////////////////////////////////////////////////////

SaveSourceData::~SaveSourceData(){
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::setup()
{
//   DataProcessingCom::setup(); // first call setup of parent class
//   DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
//   
//  
//   const CFuint nbstates = states.size();
//   sources.resize(nbstates);
//   sources_mean.resize(nbstates);  
//   
//   DIM = PhysicalModelStack::getActive()->getDim();
//   
//   for (CFuint iState = 0; iState < nbstates; ++iState)
//   {
// 	sources[iState].resize(DIM*DIM); // (source1,source2,source1_avg,source2_avg,x,y,u',v')
// 	sources_mean[iState].resize(DIM*DIM);
//   }
//   
//   for (CFuint iState = 0; iState < nbstates; ++iState)
//   {
//     CFuint nbStresses = sources[iState].size();
//     for (CFuint i = 0; i < nbStresses; ++i) {
//       sources[iState][i]=0;
//       sources_mean[iState][i]=0;
//     }
//   }
//   
//   if(m_writecoordinates) {
//     
//     std::ostringstream filenameWC;
//     filenameWC << "coordinates.dat";
//     boost::filesystem::path fnameWC = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filenameWC.str());  
//   
//     Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleWC = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();  
//   
//     ofstream& fileWC = fhandleWC->open(fnameWC);
// 
//     
//     if(DIM==2) {
//       
//       fileWC << "TITLE     =  \" Coordinates \" \n ";
//       fileWC << " VARIABLES = \"x\" \n \"y\" \n";
// 
//        for (CFuint iState = 0; iState < nbstates; ++iState) {
//          const CFreal x = (states[iState]->getCoordinates())[XX];
//          const CFreal y = (states[iState]->getCoordinates())[YY];
//          fileWC << x << "\t" << y << "\n";
//        }
//     }
// 
//     if(DIM==3) {
//       
//       fileWC << "TITLE     =  \" Coordinates \" \n ";
//       fileWC << " VARIABLES = \"x\" \t \"y\" \t \"z\" \n";
// 
//        for (CFuint iState = 0; iState < nbstates; ++iState) {
//          const CFreal x = (states[iState]->getCoordinates())[XX];
//          const CFreal y = (states[iState]->getCoordinates())[YY];
// 	 const CFreal z = (states[iState]->getCoordinates())[ZZ];
//          fileWC << x << "\t" << y << "\n";
//        }
//     }
//     
//     
//      fhandleWC->close();
//   }
// 
//   m_nameOutputFile = "sources.dat";
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SaveSourceData::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_volumes);
  result.push_back(&socket_normals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::execute(){
//
//   Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();
//   const CFuint iter = ssys_status->getNbIter();
// 
//    if (iter == 0)
//      meanCount = 0;
//    else
//      meanCount += 1;
//   
// 
//   computeReynoldsGrad();
//   prepareOutputFile();
//   
//   if(!(iter % m_saveratemean))
//      writeMean();
//   
//   
//   DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();  
//   const CFuint nbstates = states.size();
//   
//   cf_assert(m_fileRMS->isopen());
//   ofstream& fout = m_file->get();
// 
//   for (CFuint iState = 0; iState < nbstates; ++iState) {
//        RealVector currSource =sources[iState];
//        CFuint nbStresses = currSource.size();
//        for (CFuint iStress = 0; iStress < nbStresses; ++iStress) {
// 	  fout << currSource[iStress] << "\t";
//        }
//        fout << "\n";
//    }
//    
//    m_file->close();
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::unsetup(){
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::configure ( Config::ConfigArgs& args ){
  DataProcessingCom::configure( args );
/*  
  
  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);

  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  if(DIM==2) {
    /// physical model (in conservative variables)
    Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> m_varSet;
    m_varSet->setup();
        
    std::string varSetName = "Euler2DCons";
    m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());
  }
  else {
    /// physical model (in conservative variables)
    Common::SelfRegistPtr<Physics::NavierStokes::Euler3DCons> m_varSet;
    m_varSet->setup();
    std::string varSetName = "Euler3DCons";
    m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());
  }

  cf_assert(m_varSet.isNotNull());
  */

}
//////////////////////////////////////////////////////////////////////////////
void SaveSourceData::writeMean()
{
//     std::ostringstream filenameWM;
//     filenameWM << "source_mean.dat";
//     boost::filesystem::path fnameWM = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filenameWM.str());  
//   
//     Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleWM = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();  
//   
//     ofstream& fileWM = fhandleWM->open(fnameWM);
//   
//     fileWM << "TITLE     =  \" Source means \" \n ";
// 
//     CFuint nbstates=sources_mean.size();
//      
//     for (CFuint iState = 0; iState < nbstates; ++iState) {
//        RealVector currSource =sources_mean[iState];
//        CFuint nbStresses = currSource.size();
//        for (CFuint iStress = 0; iStress < nbStresses; ++iStress) {
// 	  fileWM << currSource[iStress] << "\t";
//        }
//        fileWM << "\n";
//      }
//      
//      fhandleWM->close();
//   
}

//////////////////////////////////////////////////////////////////////////////
void SaveSourceData::computeReynoldsGrad()
{
//   DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
//   DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
//   DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();
//   
//   const CFuint nbstates = states.size();
//   
//     // prepare looping over the cells
//     // find the inner trs
//   vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
//   SafePtr<TopologicalRegionSet> innerTrs=0;
//   for (CFuint i = 0; i < trs.size(); ++i) {
//     if (trs[i]->hasTag("inner")) {
//       innerTrs = trs[i];
//       break;
//     }
//   }
//   if (innerTrs==0) Common::NoSuchValueException(FromHere(),"Trs with tag 'inner' not found.");
//   
//     // prepares to loop over cells by getting the GeometricEntityPool
//   Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
//     geoBuilder = getMethodData().getStdTrsGeoBuilder();
//   StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
//   geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
// 
//   const CFuint innerNbGeos=innerTrs->getLocalNbGeoEnts();
//   
//   
//   
//   for (CFuint iState = 0; iState < nbstates; ++iState) {
//  
// // compute the gradients and the sources over the cells
//     for (CFuint iCell=0; iCell<innerNbGeos; iCell++) {
//     
// //       build the GeometricEntity
//       const CFuint cellID=innerTrs->getLocalGeoID(iCell);
//       geoData.idx = cellID;
//       GeometricEntity& cell = *geoBuilder->buildGE();
//       vector<State*> *const statesInCell = cell.getStates();
//       const CFuint nbStatesInCell = statesInCell->size();
// 
//       if(DIM==2) {
//       
//       //compute the stresses and their derivatives
//       CFreal uu = 0.0;
//       CFreal vv = 0.0;
//    
//       CFreal uv = 0.0;
//    
//       CFreal duudx = 0.0;
//       CFreal dvvdy = 0.0;
//    
//       CFreal duvdx = 0.0;
//       CFreal duvdy = 0.0;
//    
//       for (CFuint elemState = 0; elemState <nbStatesInCell; ++elemState) {
//         State *const currState = (*statesInCell)[elemState];
// 	CFuint stateID = currState->getLocalID();
//         CFreal rho = (*currState)[0];
//         CFreal u = (*currState)[1]/rho;
//         CFreal v = (*currState)[2]/rho;
//   
// 	uu = u*u;
// 	vv = v*v;
// 
// 	uv = u*v;
// 
//        	duudx +=m_normals[cellID]->getNodalNormComp(elemState,0)*uu;
// 	dvvdy +=m_normals[cellID]->getNodalNormComp(elemState,1)*vv;
// 	duvdx +=m_normals[cellID]->getNodalNormComp(elemState,0)*uv;
// 	duvdy +=m_normals[cellID]->getNodalNormComp(elemState,1)*uv;
// 
//         }
//      
//      CFreal oneover2vol = 1./(2.0*volumes[cellID]);
//    
//      duudx *= oneover2vol;
//      dvvdy *= oneover2vol;
//      duvdx *= oneover2vol;
//      duvdy *= oneover2vol;
//     
//      CFreal rho_0 = 1.0; //this is now hardcooded, need to be adopted to meanflow
//      sources[iState][0] = rho_0*duudx;
//      sources[iState][1] = rho_0*duvdy;
//      sources[iState][2] = rho_0*duvdx;
//      sources[iState][3] = rho_0*dvvdy;
//      }
//      
//      else {
//        //compute the stresses and their derivatives
//       CFreal uu = 0.0;
//       CFreal vv = 0.0;
//       CFreal ww = 0.0;
//    
//       CFreal uv = 0.0;
//       CFreal uw = 0.0;
//       CFreal vw = 0.0;
//    
//       CFreal duudx = 0.0;
//       CFreal dvvdy = 0.0;
//       CFreal dwwdz = 0.0;
//    
//       CFreal duvdx = 0.0;
//       CFreal duvdy = 0.0;
//       CFreal duwdx = 0.0;
//       CFreal duwdz = 0.0;
//       CFreal dvwdz = 0.0;
//       CFreal dvwdy = 0.0;
//       
//    
//       for (CFuint elemState = 0; elemState <nbStatesInCell; ++elemState) {
//         State *const currState = (*statesInCell)[elemState];
// 	CFuint stateID = currState->getLocalID();
//         CFreal rho = (*currState)[0];
//         CFreal u = (*currState)[1]/rho;
//         CFreal v = (*currState)[2]/rho;
//         CFreal w = (*currState)[3]/rho;
//   
// 	uu = u*u;
// 	vv = v*v;
// 	ww = w*w;
// 
// 	uv = u*v;
// 	uw = u*w;
// 	vw = v*w;
// 
//        	duudx +=m_normals[cellID]->getNodalNormComp(elemState,0)*uu;
// 	dvvdy +=m_normals[cellID]->getNodalNormComp(elemState,1)*vv;
// 	dwwdz +=m_normals[cellID]->getNodalNormComp(elemState,2)*ww;
// 	duvdx +=m_normals[cellID]->getNodalNormComp(elemState,0)*uv;
// 	duvdy +=m_normals[cellID]->getNodalNormComp(elemState,1)*uv;
// 	duwdx +=m_normals[cellID]->getNodalNormComp(elemState,0)*uw;
// 	duwdz +=m_normals[cellID]->getNodalNormComp(elemState,2)*uw;
// 	dvwdy +=m_normals[cellID]->getNodalNormComp(elemState,1)*vw;
// 	duwdz +=m_normals[cellID]->getNodalNormComp(elemState,2)*vw;
//         }
//      
//      CFreal oneover2vol = 1./(2.0*volumes[cellID]);
//    
//      duudx *= oneover2vol;
//      dvvdy *= oneover2vol;
//      dwwdz *= oneover2vol;
//      duvdx *= oneover2vol;
//      duvdy *= oneover2vol;
//      duwdx *= oneover2vol;
//      duwdz *= oneover2vol;
//      dvwdy *= oneover2vol;
//      dvwdz *= oneover2vol;
//     
//      CFreal rho_0 = 1.0; //this is now hardcooded, need to be adopted to meanflow
//      sources[iState][0] += rho_0*duudx/4.;
//      sources[iState][1] += rho_0*duvdy/4.;
//      sources[iState][2] += rho_0*duwdz/4.;
//      sources[iState][3] += rho_0*duvdx/4.;
//      sources[iState][4] += rho_0*dvvdy/4.;      
//      sources[iState][5] += rho_0*dvwdz/4.;
//      sources[iState][6] += rho_0*duwdx/4.;
//      sources[iState][7] += rho_0*dvwdy/4.;
//      sources[iState][8] += rho_0*dwwdz/4.;
//      }
//     }
//   }
  


}
//////////////////////////////////////////////////////////////////////////////
void SaveSourceData::prepareOutputFile()
{

  using boost::filesystem::path;

  cf_assert (!m_file->isopen());
  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFile);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );


  ofstream& fout = m_file->open(file);

  if(DIM==2) {
    fout << "TITLE  =  Reynold Stress Gradients in the whole field" << "\n";
    fout << "VARIABLES = dRuuDx dRuvDx dRuvDy dRvvDy " << "\n";
  }
  else {
    fout << "TITLE  =  Reynold Stress Gradients in the whole field" << "\n";
    fout << "VARIABLES = dRuuDx dRuvDx dRuwDx dRuvDy dRvvDy dRvwDy dRuwDz dRvwDz dRwwDz " << "\n";
  }

  
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////