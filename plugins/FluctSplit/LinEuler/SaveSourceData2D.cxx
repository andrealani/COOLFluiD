#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SaveSourceData2D.hh"
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

MethodCommandProvider<SaveSourceData2D,
                      DataProcessingData,
                      LinearizedEulerModule>
aSaveSourceData2DProvider("SaveSourceData2D");
//////////////////////////////////////////////////////////////////////////////

void SaveSourceData2D::defineConfigOptions(Config::OptionList& options)
{

}

//////////////////////////////////////////////////////////////////////////////

SaveSourceData2D::SaveSourceData2D( const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_sourceSave("sourceSave")
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

SaveSourceData2D::~SaveSourceData2D(){
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData2D::setup()
{
  DataProcessingCom::setup(); // first call setup of parent class
}// setup()

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SaveSourceData2D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_sourceSave);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData2D::execute(){

  //open file
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> sources = socket_sourceSave.getDataHandle();
  CFreal time =SubSystemStatusStack::getActive()->getNbIter();
  std::ostringstream filename;
  filename << "sources_"<<(time)<<".plt";
  cout << "Save sources to " << filename.str() <<".\n";
  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  Common::SelfRegistPtr<FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<FileHandlerOutput>::getInstance().create();
  ofstream& file = fhandle->open(fname);

  //write source data
  file << "TITLE = \"SourceData\"\n";
  file << "VARIABLES = \"x0\", \"x1\", \"source_x_avg\", \"source_y_avg\", \"source_x\", \"source_y\" \n";
  file << "ZONE T=\"WholeDomain\"\n";
  
  const CFuint nbcells = MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();
  
  for (CFuint iState = 0; iState < nbcells; ++iState) {
        	RealVector& sources_state = sources[iState];
		file << sources_state[0] << "\t" << sources_state[1] << "\t" <<  sources_state[2] << "\t" << sources_state[3] << "\t" <<  sources_state[4] << "\t" << sources_state[5] <<"\n";
   }
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData2D::unsetup(){
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void SaveSourceData2D::configure ( Config::ConfigArgs& args ){
  DataProcessingCom::configure( args );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////