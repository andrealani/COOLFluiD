#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SaveSourceData.hh"
#include "LinEuler/LinearizedEuler.hh"

//#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Physics {

      namespace LinearizedEuler {


//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SaveSourceData,
                      DataProcessingData,
                      LinearizedEulerModule>
aSaveSourceDataProvider("SaveSourceData");
//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::defineConfigOptions(Config::OptionList& options)
{

}

//////////////////////////////////////////////////////////////////////////////

SaveSourceData::SaveSourceData( const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_sources("sources")
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

SaveSourceData::~SaveSourceData(){
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::setup()
{
  DataProcessingCom::setup(); // first call setup of parent class
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> sources = socket_sources.getDataHandle();
  sources.resize(states.size());
  for (CFuint iState = 0; iState < sources.size(); ++iState)
  {
	sources[iState].resize(8); // (source1,source2,source1_avg,source2_avg,x,y,u',v')
  }
  for (CFuint iState = 0; iState < sources.size(); ++iState)
  {
	sources[iState]=0;
  }
}// setup()

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SaveSourceData::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
SaveSourceData::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_sources);
  return result;
}
//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::execute(){

  CFAUTOTRACE;

  //open file
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> sources = socket_sources.getDataHandle();
  CFreal time =SubSystemStatusStack::getActive()->getNbIter();
  std::ostringstream filename;
  filename << "sources_"<<(time)<<".plt";
  cout << "Save sources to " << filename.str() <<".\n";
  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  Common::SelfRegistPtr<FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<FileHandlerOutput>::getInstance().create();
  ofstream& file = fhandle->open(fname);

  //write source data
  file << "TITLE = \"SourceData\"\n";
  file << "VARIABLES = \"x0\", \"x1\", \"source_x\", \"source_y\", \"source_x_avg\", \"source_y_avg\", \"u\", \"v\"\n";
  file << "ZONE T=\"WholeDomain\"\n";
  
  const CFuint nbcells = MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();
  for (CFuint iState = 0; iState < nbcells; ++iState) {
        	RealVector& sources_state = sources[iState];
		file << sources_state[4] << "\t" << sources_state[5] << "\t" <<  sources_state[0] << "\t" << sources_state[1] << "\t" <<  sources_state[2] << "\t" << sources_state[3] << "\t"<< sources_state[6] << "\t" << sources_state[7] <<"\n";
   }
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::unsetup(){
  DataHandle<RealVector> sources = socket_sources.getDataHandle();
  sources.resize(0);
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData::configure ( Config::ConfigArgs& args ){
  DataProcessingCom::configure( args );

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////