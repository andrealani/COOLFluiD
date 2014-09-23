#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"

#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/StatesPeriodicBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace std;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StatesPeriodicBC,FiniteElementMethodData,FiniteElementModule >
  statesPeriodicBCProvider("StatesPeriodicBC");

//////////////////////////////////////////////////////////////////////////////

void StatesPeriodicBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >( "FileName",
     "Name of file containing list of states to apply dirichlet BC on." );
   options.addConfigOption< std::string >( "OtherFileName",
     "Name of file containing list of states to apply dirichlet BC on." );
   options.addConfigOption< CFreal >( "RotationAngle",
     "Rotation around X-axis (in degrees)." );

   options.addConfigOption< std::vector< CFuint > >( "ApplyEqs",
     "Apply the BC only to specified equations zero-based indexed (default all equations)." );
   
}

//////////////////////////////////////////////////////////////////////////////

StatesPeriodicBC::StatesPeriodicBC(const std::string& name) :
  FiniteElementMethodCom(name),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_states("states"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_appliedStrongBC("appliedStrongBC")
{
  addConfigOptionsTo(this);

  m_functions = std::vector< std::string >();
  m_vars = std::vector< std::string >();
  m_applyEqs = std::vector< CFuint >();

  setParameter( "ApplyEqs", &m_applyEqs);

  setParameter( "FileName", &m_filename);
  setParameter( "OtherFileName", &m_otherfilename);
  setParameter( "RotationAngle", &m_rotationAngle);

}

//////////////////////////////////////////////////////////////////////////////

void StatesPeriodicBC::setup()
{
  CFAUTOTRACE;

  // validate ApplyEqs config option
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (!m_applyEqs.size()) {
    // apply to all PhysicalModel equations
    m_applyEqs.resize(nbEqs);
    for (CFuint i=0; i<nbEqs; ++i)
      m_applyEqs[i]=i;
  } else {
    // validate the equations to apply the BC exist
    for (CFuint i=0; i<m_applyEqs.size(); ++i)
      if (m_applyEqs[i]>= nbEqs)
        throw BadValueException(
          "StatesPeriodicBC ApplyEqs refers to an equation that doesn't exist." );
  }

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

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
  CFout << "Finished reading boundary file\n";

  fhandle->close();


///Read File and build list of States
  using namespace boost::filesystem;

  // Read the file and fill in the statesList vector
  path file2 = Environment::DirPaths::getInstance().getWorkingDir() / path(m_otherfilename);

  SelfRegistPtr<Environment::FileHandlerInput> fhandle2 = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin2 = fhandle2->open(file2);

  //Check content of first line
  getline(fin2,line);
  words = Common::StringOps::getWords(line);
  if(words.size() != 1) {
    throw BadFormatException("Wrong number of parameters in 1st line of file: " + file.string());
  }

  // Check number of nodes identifier
  if(words[0] != "!COOLFLUID_DIRICHLET"){
    throw BadFormatException("Expecting !COOLFLUID_DIRICHLET identifier in 1st line of " + file.string());
  }

  //Check content of second line
  getline(fin2,line);
  words = Common::StringOps::getWords(line);
  if(words.size() != 2)  {
    throw BadFormatException("Wrong number of parameters in 2nd line of file: " + file.string());
  }

  // Check number of states identifier
  if(words[0] != "!NB_STATES"){
    throw BadFormatException("Expecting !NBSTATES identifier in 2nd line of " + file.string());
  }
  // Check agreement of the number of states
  CFuint nbDirichletStates2 = from_str<CFint>(words[1]);
  m_otherStatesList.resize(nbDirichletStates2);

  //Read the list of States
  for (CFuint iState = 0; iState < nbDirichletStates2; ++iState) {
    fin2 >> m_otherStatesList[iState];
  }
  CFout << "Finished reading matching boundary file\n";

  fhandle2->close();

  assert(nbDirichletStates == nbDirichletStates2);
  m_correspondingStatesList.resize(m_otherStatesList.size());

  for (CFuint iState = 0; iState < nbDirichletStates; ++iState) {
    const CFuint nLocalID = m_statesList[iState];
    assert(nLocalID < states.size());
    const State *currState = states[nLocalID];
    RealVector currNode = states[nLocalID]->getCoordinates();
    CFreal minimumDistance = MathTools::MathConsts::CFrealMax();
    bool foundMatch = false;
    CFreal threshold = 0.0001;
    RealVector coord(DIM_3D);

    for (CFuint iOtherState = 0; iOtherState < nbDirichletStates; ++iOtherState) {
      const CFuint otherLocalID = m_otherStatesList[iOtherState];
      assert(otherLocalID < states.size());
      const State *otherState = states[otherLocalID];
      RealVector otherNode = states[otherLocalID]->getCoordinates();

      ///Rotate the node around X-axis
      CFreal cosT = cos(m_rotationAngle*MathTools::MathConsts::CFrealPi()/180.);
      CFreal sinT = sin(m_rotationAngle*MathTools::MathConsts::CFrealPi()/180.);

      for (CFuint j = 0; j < DIM_3D; ++j) {
        // Get the old coordinates
        coord[j] = otherNode[j];
      }

      // Modify the coordinates: Rotation
      coord[0] = coord[0];
      coord[1] = cosT*(coord[1]) + sinT*(coord[2]);
      coord[2] = -sinT*(coord[1]) + cosT*(coord[2]);

      for (CFuint j = 0; j < DIM_2D; ++j) {
        // Set the new coordinates
        otherNode[j] = coord[j];
      }

      ///Compute distance between rotated point and node
      CFreal distance = MathTools::MathFunctions::getDistance(currNode, otherNode);

      if (distance < minimumDistance){
        minimumDistance = distance;
        m_correspondingStatesList[iState] = otherLocalID;
        if(distance < threshold) foundMatch = true;
      }
    }
    CFout << "minimumDistance: " << minimumDistance << "\n";
    assert(foundMatch == true);
  }

}

//////////////////////////////////////////////////////////////////////////////

void StatesPeriodicBC::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal > rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool > isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<std::valarray< State* > > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();

  // get system matrix and global index mapping
  const Common::SafePtr< LinearSystemSolver > lss =
    getMethodData().getLinearSystemSolver()[0];
  Common::SafePtr< LSSMatrix > sysMat = lss->getMatrix();
  const LSSIdxMapping& idxMapping = lss->getLocalToGlobalMapping();

  // PhysicalModel properties and auxiliary variables
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();
  const CFuint nbApplyEqs = m_applyEqs.size();

  // number of variables: dimensions, time and PhysicalModel dimensions
  RealVector dirichletState(nbEqs);
  RealVector dirichletRhs(nbApplyEqs);

  // this should be an intermediate lightweight assembly! it is needed because
  // here you SET values while elsewhere you ADD values
  sysMat->flushAssembly();

  // cycle all the states in the TRS
  for(CFuint iState = 0; iState < m_statesList.size(); ++iState) {
    const CFuint nLocalID = m_statesList[iState];
    assert(nLocalID < states.size());
    const State *currState = states[nLocalID];
    if (isUpdated[nLocalID]) continue;
    if (currState->isParUpdatable()) {

      // global position of node, which must have at least one neighbour
      const CFuint nGlobalID = idxMapping.getColID(nLocalID)*nbEqs;
      const CFuint nbNeigh = bStatesNeighbors[nLocalID].size();
      cf_assert(nbNeigh>0);

      for(CFuint j=0; j<nbApplyEqs; ++j) {
        const CFuint jEq=m_applyEqs[j];
        dirichletRhs[j] = 0.;
      }


      // erase system matrix line (all its columns, including the diagonal)
      // unless symmetry method is ScaleDiagonal. also, if symmetry method
      // is AdjustColumn, pass column contribution to the rhs vector then
      // erase all its lines from system matrix. only cycle on neighbour
      // nodes to avoid expensive reallocations
      for (CFuint j=0; j<nbNeigh; ++j) {
        // global position of neighbour node
        const CFuint jGlobalID = idxMapping.getColID(
          bStatesNeighbors[nLocalID][j]->getLocalID() )*nbEqs;

        // for the specific node equation line where to apply the BC,
        // erase the line (all its columns)
        for (CFuint i=0; i<nbApplyEqs; ++i) {
          const CFuint iEq=m_applyEqs[i];
          for (CFuint jEq=0; jEq<nbEqs; ++jEq)
            sysMat->setValue(nGlobalID+iEq,jGlobalID+jEq, 0.);
        }
      } // for nbNeigh

      //get the ID of the corresponding state
      const CFuint otherLocalID = m_correspondingStatesList[iState];
      const CFuint otherGlobalID = idxMapping.getColID(otherLocalID)*nbEqs;

      // set rhs and system matrix diagonal terms (scaled)
      for (CFuint j=0; j<nbApplyEqs; ++j) {
        const CFuint jEq=m_applyEqs[j];
        sysMat->setValue(nGlobalID+jEq,otherGlobalID+jEq, -1.);
        sysMat->setValue(nGlobalID+jEq,nGlobalID+jEq, 1.);
        rhs(nLocalID,jEq,nbEqs) = dirichletRhs[j];
      }

      // flagging is important!
      isUpdated[nLocalID] = true;

    } // isParUpdatable?

  } // cycle all the states in the StatesList

  // this should be an intermediate lightweight assembly! it is needed because
  // here you SET values while elsewhere you ADD values
  sysMat->flushAssembly();

}

//////////////////////////////////////////////////////////////////////////////

void StatesPeriodicBC::configure ( Config::ConfigArgs& args )
{
  FiniteElementMethodCom::configure(args);

  m_vFunction.setFunctions(m_functions);
  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StatesPeriodicBC::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_states);
  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_appliedStrongBC);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

