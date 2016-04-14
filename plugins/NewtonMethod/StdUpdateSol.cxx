// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "StdUpdateSol.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/ConvergenceStatus.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonIterator.hh"

#include "Framework/PathAppender.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUpdateSol, NewtonIteratorData, NewtonMethodModule> 
stdUpdateSolProvider("StdUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSol::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<vector<CFreal>,Config::DynamicOption<> >("Relaxation","Relaxation factor");
  options.addConfigOption< bool >("Validate","Check that each update creates variables with physical meaning");
}

//////////////////////////////////////////////////////////////////////////////

StdUpdateSol::StdUpdateSol(const std::string& name) : 
  NewtonIteratorCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_invalidStates("invalidStates")
{
  addConfigOptionsTo(this);
   
  m_alpha = vector<CFreal>();
  setParameter("Relaxation",&m_alpha);

  m_validate = false;
  setParameter("Validate",&m_validate);
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSol::setup()
{
  const CFuint nbStates = socket_states.getDataHandle().size();

  if (m_validate) {
    socket_invalidStates.getDataHandle().resize(nbStates);
  }
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (m_alpha.size() == 0) {
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = 1.;
    }
  }

  if (m_alpha.size() == 1) {
    const CFreal value = m_alpha[0];
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = value;
    }
  }

  if (m_alpha.size() != nbEqs) {
    throw BadValueException (FromHere(),"StdUpdateSol::setup() : m_alpha.size() != nbEqs");
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSol::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint states_size = states.size();

  SafePtr<FilterState> filterState = getMethodData().getFilterState();
  SafePtr<FilterRHS> filterRHS     = getMethodData().getFilterRHS();
  
  for (CFuint iState = 0; iState < states_size; ++iState)
  {
    State& cur_state = *states[iState];
    // do the update only if the state is parallel updatable
    if (cur_state.isParUpdatable())
    {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	filterRHS->filter(iEq, dU(iState, iEq, nbEqs));
	cur_state[iEq] += m_alpha[iEq] * dU(iState, iEq, nbEqs);
      }
      
      // apply a polymorphic filter to the state
      filterState->filter(cur_state);
    }
    else {
      // reset to 0 the RHS for ghost states in order to avoid 
      // inconsistencies in the parallel L2 norm computation
      cf_assert(!cur_state.isParUpdatable());
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	dU(iState, iEq, nbEqs) = 0.;
      }
    }
  }
  
  ///This vector will temporally contain the ID of any vector considered not valid.
  std::vector<CFuint> badStatesIDs;

  ///Name of the file where to store invalid states in.
  boost::filesystem::path filepath ( "unphysical_states.plt");

  filepath = PathAppender::getInstance().appendParallel( filepath );

  // loop for unphysicalness check:
  if ( m_validate )
  {
    Common::SafePtr<SpaceMethod> theSpaceMethod = getMethodData().getCollaborator<SpaceMethod>();
    Common::SafePtr<SpaceMethodData> theSpaceMethodData = theSpaceMethod->getSpaceMethodData();
    Common::SafePtr<Framework::ConvectiveVarSet> theVarSet = theSpaceMethodData->getUpdateVar();

    ///Reset the vector:
    badStatesIDs.resize(0);

    for (CFuint iState = 0; iState < states_size; ++iState)
    {
      State& cur_state = *states[iState];
      if ( !( theVarSet->isValid(cur_state) ) )
      {
        badStatesIDs.push_back( iState );
      }
    } // for all states
  } // if validate

  if (badStatesIDs.size() > 0)
  {
    cf_assert(m_validate);
    // sort the invalid states in order to be able to search later on 
    sort(badStatesIDs.begin(), badStatesIDs.end());
    
    DataHandle<CFreal> invalidStates = socket_invalidStates.getDataHandle();
    cf_assert(invalidStates.size() == states_size);
    invalidStates = 0.;
    for (CFuint i = 0; i < badStatesIDs.size(); ++i) {
      invalidStates[badStatesIDs[i]] = 1.;
    }
    
    CFLog(VERBOSE, "StdUpdateSol::execute() => [" << badStatesIDs.size() << "] invalid states detected\n");
    // correct all unphysical states 
    correctUnphysicalStates(badStatesIDs);
    
    ///Gets the global iteration
    CFuint cur_global_iter = SubSystemStatusStack::getActive()->getNbIter();
    ///Gets the Newton method procedure iteration
    CFuint cur_iter = getMethodData().getConvergenceStatus().iter;

    ///opens a file to put the states which have unphysical variables
    ofstream fout ( filepath.string().c_str(), ios::app );

    if ( !fout.is_open() )
    {
      cout << endl << "File [" << filepath.string().c_str() << "]  is not open.\n";
    }
    ///
    const std::vector<std::string>& theVarNames=getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVar()->getVarNames();
    CFuint dim = PhysicalModelStack::getActive()->getDim();
    ///

    ///Zone header:
    fout << "VARIABLES = \"X\"\n ";
    fout << "\"Y\"" << endl;
    if (dim == 3 ){
        fout << "\"Z\"" << endl;
    }
    for (CFuint i = 0; i < theVarNames.size(); i++){
      fout << "\"" << theVarNames[i] << "\"" << endl;
    }
    ///


    CFuint badStates_size = badStatesIDs.size();

    ///Line  
    fout << "ZONE T = \"Global iteration "<< cur_global_iter << ", Newton iteration " << cur_iter << "\"" << endl;
    ///Line
    fout << "I="<< badStates_size << ", J=1, K=1, ZONETYPE=Ordered"<< endl;
    ///Line
    fout << "DATAPACKING=POINT" << endl;
    ///Line
    fout << "DT=(";
    for ( CFuint i = 0 ; i < dim + theVarNames.size(); i++){
      fout << "DOUBLE ";
    }

    fout << ")" << endl;
    ///

    for (CFuint index = 0; index < badStates_size; ++index)
    {
     State& invalid_state  = *states[ badStatesIDs[index] ];
     fout << invalid_state.getCoordinates() << " " << invalid_state << std::endl;
    } // for all invalid states

    fout.close();

  } // writing badStates

  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdUpdateSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> > StdUpdateSol::providesSockets()
{
  vector<SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_invalidStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
