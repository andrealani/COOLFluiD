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
#include "FiniteElement/StatesDirichletBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace std;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StatesDirichletBC,FiniteElementMethodData,FiniteElementModule >
  statesDirichletBCProvider("StatesDirichletBC");

//////////////////////////////////////////////////////////////////////////////

void StatesDirichletBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >( "FileName",
     "Name of file containing list of states to apply dirichlet BC on." );


   options.addConfigOption< std::string >( "Symmetry",
     "Keep matrix symmetry by AdjustColumn or ScaleDiagonal methods or not (default)." );
   options.addConfigOption< std::vector< std::string > >( "Vars",
     "Definition of the Variables (required)." );
   options.addConfigOption< std::vector< CFuint > >( "ApplyEqs",
     "Apply the BC only to specified equations zero-based indexed (default all equations)." );
   options.addConfigOption< CFreal >( "ScaleDiagonal",
     "ScaleDiagonal symmetry method coefficient (default 1.e20)." );
   options.addConfigOption< std::vector< std::string > >( "Def",
     "Definition of the Functions (required)." );
   options.addConfigOption< bool >( "Implicit",
     "Apply the BC implicitly? (default false)" );
}

//////////////////////////////////////////////////////////////////////////////

StatesDirichletBC::StatesDirichletBC(const std::string& name) :
  FiniteElementMethodCom(name),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_states("states"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_appliedStrongBC("appliedStrongBC")
{
  addConfigOptionsTo(this);

  m_isImplicit = false;
  m_functions = std::vector< std::string >();
  m_vars = std::vector< std::string >();
  m_applyEqs = std::vector< CFuint >();
  m_symmetryStr = "No";
  m_scale = 1.e20;

  setParameter( "Implicit",      &m_isImplicit);
  setParameter( "Def",           &m_functions);
  setParameter( "Vars",          &m_vars);
  setParameter( "ApplyEqs",      &m_applyEqs);
  setParameter( "Symmetry",      &m_symmetryStr);
  setParameter( "ScaleDiagonal", &m_scale);

  setParameter( "FileName",      &m_filename);

}

//////////////////////////////////////////////////////////////////////////////

void StatesDirichletBC::setup()
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
          "StatesDirichletBC ApplyEqs refers to an equation that doesn't exist." );
  }


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

  ///prepare the isAppliedStrongBC datahandle
  DataHandle< vector<bool> > appliedStrongBC =
    socket_appliedStrongBC.getDataHandle();

  for(CFuint iState = 0; iState < m_statesList.size(); ++iState) {
        const CFuint localStateID = m_statesList[iState];
	
	appliedStrongBC[localStateID].resize(nbEqs);
	for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
		appliedStrongBC[localStateID][iVar] =
		(m_applyEqs[iVar]) ? true : false;
	}
  }




}

//////////////////////////////////////////////////////////////////////////////

void StatesDirichletBC::execute()
{
  CFAUTOTRACE;

  getMethodData().setDirichletBCApplied(true);

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
  const bool isAdjust = (m_symmetryStr=="AdjustColumn");
  const bool isScale  = (m_symmetryStr=="ScaleDiagonal");
  if (!isScale)
    m_scale = 1.;
  CFreal implicit = (m_isImplicit? 1.:0.);


  // number of variables: dimensions, time and PhysicalModel dimensions
  RealVector dirichletState(nbEqs);
  RealVector dirichletRhs(nbApplyEqs);
  RealVector variables(dim+1+nbEqs);
  variables[dim] = SubSystemStatusStack::getActive()->getCurrentTimeDim();


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

      // evaluate boundary condition at space, time (omitted) and states,
      // and calculate the boundary condition enforced value
      const RealVector& temp = currState->getCoordinates();
      for (CFuint i=0; i<temp.size(); ++i)
        variables[i] = temp[i];
      for(CFuint i=0; i<nbEqs; ++i)
        variables[dim+1+i] = (*currState)[i];
      m_vFunction.evaluate(variables,dirichletState);
// CFout << "DirichletState Assigning value: " << dirichletState <<  " to state: " << currState->getLocalID() << " at coord: " << currState->getCoordinates() << "\n";

      for(CFuint j=0; j<nbApplyEqs; ++j) {
        const CFuint jEq=m_applyEqs[j];
        dirichletRhs[j] =
          dirichletState[jEq] - implicit*((*currState)[jEq]);
      }


      // erase system matrix line (all its columns, including the diagonal)
      // unless symmetry method is ScaleDiagonal. also, if symmetry method
      // is AdjustColumn, pass column contribution to the rhs vector then
      // erase all its lines from system matrix. only cycle on neighbour
      // nodes to avoid expensive reallocations
      if (!isScale) {

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

        if (isAdjust) {

          // set affecting rows and the columns ids
          const CFuint M = nbApplyEqs*nbNeigh;
          const CFuint N = nbApplyEqs;
          CFint iM[M], iN[N];
          for (CFuint i=0; i<nbNeigh; ++i) {
            const CFuint iGlobalID = idxMapping.getColID(
              bStatesNeighbors[nLocalID][i]->getLocalID() )*nbEqs;
            for (CFuint j=0; j<nbApplyEqs; ++j)
              iM[i*nbApplyEqs+j] = iGlobalID + m_applyEqs[j];
          }
          for (CFuint j=0; j<nbApplyEqs; ++j)
            iN[j] = nGlobalID + m_applyEqs[j];

          // get matrix terms with node's neighbour's contributions and
          // transfer their contribution to rhs
          RealMatrix sysMatM(M,N);
          sysMat->finalAssembly();
          sysMat->getValues( M,iM, N,iN, &sysMatM[0] );
          for (CFuint i=0; i<nbNeigh; ++i) {
            const CFuint neighStateID =
              bStatesNeighbors[nLocalID][i]->getLocalID();
            for (CFuint j=0; j<nbApplyEqs; ++j) {
              const CFuint jEq = m_applyEqs[j];
              for (CFuint k=0; k<nbApplyEqs; ++k)
                rhs(neighStateID,jEq,nbEqs) -=
                  sysMatM(i*nbApplyEqs+j,k) * dirichletRhs[k];
            }
          }
          sysMatM = 0.;
          sysMat->setValues( M,iM, N,iN, &sysMatM[0] );

        } // isAdjust?
      } // !isScale?


      // set rhs and system matrix diagonal terms (scaled)
      for (CFuint j=0; j<nbApplyEqs; ++j) {
        const CFuint jEq=m_applyEqs[j];
        sysMat->setValue(nGlobalID+jEq,nGlobalID+jEq, m_scale);
        rhs(nLocalID,jEq,nbEqs) = dirichletRhs[j] * m_scale;
      }

      // flagging is important!
      isUpdated[nLocalID] = true;

    } // isParUpdatable?

  } // cycle all the states in the TRS

  // this should be an intermediate lightweight assembly! it is needed because
  // here you SET values while elsewhere you ADD values
  sysMat->flushAssembly();

}

//////////////////////////////////////////////////////////////////////////////

void StatesDirichletBC::configure ( Config::ConfigArgs& args )
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
  StatesDirichletBC::needsSockets()
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

