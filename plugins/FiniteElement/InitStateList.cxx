#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"

#include "Common/CFLog.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/InitStateList.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateList, FiniteElementMethodData, FiniteElementModule> initStateListProvider("InitStateList");

//////////////////////////////////////////////////////////////////////////////

void InitStateList::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >( "FileName",
     "Name of file containing list of states to apply dirichlet BC on." );

   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

InitStateList::InitStateList(const std::string& name) :
  FiniteElementMethodCom(name),
  socket_states("states")
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);

  setParameter( "FileName",      &m_filename);

}

//////////////////////////////////////////////////////////////////////////////

InitStateList::~InitStateList()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitStateList::setup()
{
  CFAUTOTRACE;

///Read File and build list of States
  using namespace boost::filesystem;

  // Read the file and fill in the statesList vector
  path file = Environment::DirPaths::getInstance().getWorkingDir() / path(m_filename);

  SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(file);

  std::string line;
  vector<std::string> words;

  //Check content of first line
  getline(fin,line);
  words = Common::StringOps::getWords(line);
  if(words.size() != 1) {
    throw BadFormatException("Wrong number of parameters in 1st line of file: " + file.string());
  }

  // Check number of nodes identifier
  if(words[0] != "!COOLFLUID_DIRICHLET"){
    throw BadFormatException("Expecting !COOLFLUID_DIRICHLET identifier in 1st line of " + file.string());
  }

  //Check content of second line
  getline(fin,line);
  words = Common::StringOps::getWords(line);
  if(words.size() != 2)  {
    throw BadFormatException("Wrong number of parameters in 2nd line of file: " + file.string());
  }

  // Check number of states identifier
  if(words[0] != "!NB_STATES"){
    throw BadFormatException("Expecting !NBSTATES identifier in 2nd line of " + file.string());
  }
  // Check agreement of the number of states
  CFuint nbDirichletStates = from_str<CFint>(words[1]);
  m_statesList.resize(nbDirichletStates);

  //Read the list of States
  for (CFuint iState = 0; iState < nbDirichletStates; ++iState) {
    fin >> m_statesList[iState];
  }
CFout << "Read statesID finished\n";

  fhandle->close();

}

//////////////////////////////////////////////////////////////////////////////

void InitStateList::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FiniteElementMethodCom::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void InitStateList::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealVector dirichletState(nbEqs);

  for(CFuint iState = 0; iState < m_statesList.size(); ++iState) {
    const CFuint nLocalID = m_statesList[iState];
    assert(nLocalID < states.size());
    State *currState = states[nLocalID];

    const RealVector& node = currState->getCoordinates();
      _vFunction.evaluate(node,dirichletState);
      *currState = dirichletState;
//     CFout << "Assigning state: " << states[*itd]->getLocalID() << " at coord: " << states[*itd]->getCoordinates() << " the value: " << *(states[*itd]) <<"\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
InitStateList::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
