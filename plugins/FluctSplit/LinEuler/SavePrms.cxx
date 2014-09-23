#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SavePrms.hh"
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

MethodCommandProvider<SavePrms,
                      DataProcessingData,
                      LinearizedEulerModule>
aSavePrmsProvider("SavePrms");

//////////////////////////////////////////////////////////////////////////////

SavePrms::SavePrms( const std::string& name) :
  DataProcessingCom(name),
  socket_states("states")
{
  addConfigOptionsTo(this);
  
  m_saveRate = 1;
  setParameter("SaveRate",&m_saveRate);
}

//////////////////////////////////////////////////////////////////////////////

SavePrms::~SavePrms(){
}
//////////////////////////////////////////////////////////////////////////////

void SavePrms::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file with."); 
}

//////////////////////////////////////////////////////////////////////////////

void SavePrms::setup()
{
  DataProcessingCom::setup(); // first call setup of parent class
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbstates = states.size();
  
  Prms.resize(nbstates);
  for(int i=0; i<nbstates; i++)
    Prms[i] = 0.;
}// setup()

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SavePrms::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SavePrms::execute(){

  //open file
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  CFuint time =SubSystemStatusStack::getActive()->getNbIter();

  const CFuint nbstates = states.size();
  
  for (CFuint iState = 0; iState < nbstates; ++iState) {
    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    State& curr_state = *states[iState];
    const CFreal p2 = curr_state[3]*curr_state[3];
    Prms[iState] += p2;
   }

  if(!(time % m_saveRate)) {
   
    
  std::ostringstream filename;
  filename << "prms_"<<(time)<<".plt";
  cout << "Save prms to " << filename.str() <<".\n";
  boost::filesystem::path fname = PathAppender::getInstance().appendParallel(filename.str());
  Common::SelfRegistPtr<FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<FileHandlerOutput>::getInstance().create();

  
  ofstream& file = fhandle->open(fname);

  
  //write source data
  file << "TITLE = \"SourceData\"\n";
  file << "VARIABLES = \"x0\", \"x1\", \"Prms\" \n";
  file << "ZONE T=\"WholeDomain\"\n";
  
  const CFuint nbstates = states.size();
  
  for (CFuint iState = 0; iState < nbstates; ++iState) {
    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    file << x << "\t" << y << "\t" << Prms[iState] << "\n";
   }
   
  fhandle->close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SavePrms::unsetup(){
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void SavePrms::configure ( Config::ConfigArgs& args ){
  DataProcessingCom::configure( args );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////