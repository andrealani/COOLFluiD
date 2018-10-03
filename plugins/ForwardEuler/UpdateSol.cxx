// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "MathTools/MathChecks.hh"

#include "Framework/OutputFormatter.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/PathAppender.hh"
#include "Common/MPI/MPIStructDef.hh"

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/UpdateSol.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ForwardEuler {
  
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateSol,
                      FwdEulerData,
                      ForwardEulerLib>
aUpdateSolProvider("StdUpdateSol");

//////////////////////////////////////////////////////////////////////////////
      
void UpdateSol::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("Validate","Check that each update creates variables with physical meaning");

  options.addConfigOption< bool >
    ("ClipResidual","Clips residual to 0 if < 1e-16 (useful in saome cases)");
  
  options.addConfigOption< vector<string> >
    ("TrsWithNullRHS","Names of the TRSs whose states have a RHS to be set to 0");
}
    
//////////////////////////////////////////////////////////////////////////////
    
UpdateSol::UpdateSol(const std::string& name) :
  FwdEulerCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_volumes("volumes", false),
  m_activeStates()
{
  addConfigOptionsTo(this);
  m_validate = false;
  setParameter("Validate",&m_validate);
  
  m_clipResidual = false;
  setParameter("ClipResidual",&m_clipResidual);
  
  m_trsWithNullRHS = vector<string>();
  setParameter("TrsWithNullRHS",&m_trsWithNullRHS);
}
    
//////////////////////////////////////////////////////////////////////////////
     
void UpdateSol::execute()
{
  CFAUTOTRACE;

  SafePtr<FilterState> filterState = getMethodData().getFilterState();

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs  = socket_rhs.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFreal sys_dt  = SubSystemStatusStack::getActive()->getDT();
  const CFreal CFL = getMethodData().getCFL()->getCFLValue();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();
  const bool isTimeAccurate = getMethodData().isTimeAccurate();
  const bool isGlobalTimeStep = getMethodData().isGlobalTimeStep();
  bool isTimeStepTooLarge = false;
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  CFreal maxUpdateCoeff = 0.;
  if(isGlobalTimeStep){
    for (CFuint i = 0; i < nbStates; ++i) {
      if (states[i]->isParUpdatable()) {
	cf_assert(updateCoeff[i] > 0.);
	maxUpdateCoeff = std::max(maxUpdateCoeff, updateCoeff[i]);
      }
    }
    cf_assert(maxUpdateCoeff > 0.);
  }
  
  CFreal dt = 0.;
  CFreal maxCFL = 1.;
  for (CFuint i = 0; i < nbStates; ++i)
  {
    // cout << i << " => "; cout.precision(12); cout << updateCoeff[i] << endl;
      State& cur_state = *states[i];

    // do the update only if the state is parallel updatable
    if (cur_state.isParUpdatable())
    {
      const bool isZeroUpdateCoeff = MathChecks::isZero(updateCoeff[i]);
      CFLogDebugMax( "updateCoeff[i] = " << updateCoeff[i]
        << (isZeroUpdateCoeff? " (ZERO)\n" : "\n") );

      if (!isGlobalTimeStep) {
	if(!isTimeAccurate)  // non time accurate update
	  {
	    dt = (isZeroUpdateCoeff? 0. : CFL / updateCoeff[i]);
	    if (isZeroUpdateCoeff) {
	      for (CFuint j = 0; j < nbEqs; ++j)
            if (MathChecks::isNotZero(rhs(i, j, nbEqs))) {
              CFLog(VERBOSE,"Residual contribution to state "
                    << i << " with null characteristic speed.\n");
              break;
            }
	    }
	  }
	else // if time accurate update
	  {
	    // Compute maximum DT
	    const CFreal dtmax = 1./updateCoeff[i];
	    dt = sys_dt / volumes[i];
	    
	    // Compute equivalent CFL
	    const CFreal ratio = dt/dtmax;
	    
	    if(ratio > 1.)
	      {
		isTimeStepTooLarge = true;
		maxCFL = max(maxCFL,ratio);
	      }
	  }
      }
      else {
	// AL: this assumes same volume... needs to be generalized
	dt = CFL/maxUpdateCoeff;
      }
      
      CFLogDebugMax( "dt = " << dt << "\n");
      
      CFLogDebugMed("rhs = ");
      // update of the state values
      for (CFuint j = 0; j < nbEqs; ++j)
      {
	// AL: important fix for avoiding unpredictable behaviour due to small numbers ( < machine precision) in rhs
	//     which sometimes lead to floating exception while computing the norm 
	if (m_clipResidual && std::abs(rhs(i,j,nbEqs)) < 1e-16) {
	  rhs(i,j,nbEqs) = 0.;
	}
	CFLogDebugMed(rhs(i,j,nbEqs) << " ");
	
	CFreal rhsState = (m_activeStates[cur_state.getLocalID()]) ? rhs(i,j,nbEqs) : 0.;
	cur_state[j] += rhsState *= dt;
      }
      CFLogDebugMed("\n");
      
      filterState->filter(cur_state);
      cf_assert(cur_state.isValid());
    }
    else {
      // reset to 0 the RHS for ghost states in order to avoid 
      // inconsistencies in the parallel L2 norm computation
      cf_assert(!cur_state.isParUpdatable());
      for (CFuint j = 0; j < nbEqs; ++j) {
	rhs(i,j,nbEqs) = 0.;
      }
    }
    
    updateCoeff[i] = 0.0;
    
    //  std::cout.precision(12); std::cout << i << " => "<< cur_state <<"\n";
  }
    
  // This vector will temporally contain the ID of any vector considered not valid.
  std::vector<CFuint> badStatesIDs;
  
  string fileName = "unphysical_states.plt" + StringOps::to_str(PE::GetPE().GetRank("Default"));
  // Name of the file where to store invalid states in.
  boost::filesystem::path filepath (fileName.c_str());
  
  // loop for unphysicalness check:
  if ( m_validate )
  {
    filepath = PathAppender::getInstance().appendParallel( filepath );
    
    Common::SafePtr<SpaceMethod> theSpaceMethod = getMethodData().getCollaborator<SpaceMethod>();
    Common::SafePtr<SpaceMethodData> theSpaceMethodData = theSpaceMethod->getSpaceMethodData();
    Common::SafePtr<Framework::ConvectiveVarSet> theVarSet = theSpaceMethodData->getUpdateVar();
    
    // Reset the vector:
    badStatesIDs.resize(0);
    
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {
      State& cur_state = *states[iState];
      if ( !( theVarSet->isValid(cur_state) ) )
	{
	  badStatesIDs.push_back( iState );
	}
      
    } // for all states
  } // if validate
  
  if ( badStatesIDs.size() > 0 )
  {
    // Gets the global iteration
    CFuint cur_global_iter = SubSystemStatusStack::getActive()->getNbIter();
    // Gets the Newton method procedure iteration
    CFuint cur_iter = getMethodData().getConvergenceStatus().iter;
    
    // opens a file to put the states which have unphysical variables
    ofstream fout ( filepath.string().c_str(), ios::app );
    
    if ( !fout.is_open() )
      {
	cout << endl << "File [" << filepath.string().c_str() << "]  is not open.\n";
      }
    
      const std::vector<std::string>& theVarNames=getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVar()->getVarNames();
      CFuint dim = PhysicalModelStack::getActive()->getDim();
      
      //Zone header:
      fout << "VARIABLES = \"X\"\n ";
      fout << "\"Y\"" << endl;
      if (dim == 3 ){
	fout << "\"Z\"" << endl;
      }
      for (CFuint i = 0; i < theVarNames.size(); i++){
	fout << "\"" << theVarNames[i] << "\"" << endl;
      }
      
      CFuint badStates_size = badStatesIDs.size();
      
      //Line  
      fout << "ZONE T = \"Global iteration "<< cur_global_iter << ", Newton iteration " << cur_iter << "\"" << endl;
      //Line
      fout << "I="<< badStates_size << ", J=1, K=1, ZONETYPE=Ordered"<< endl;
      //Line
      fout << "DATAPACKING=POINT" << endl;
      //Line
      fout << "DT=(";
      for ( CFuint i = 0 ; i < dim + theVarNames.size(); i++){
	fout << "DOUBLE ";
      }
      
      fout << ")" << endl;
      //
      
      for (CFuint index = 0; index < badStates_size; ++index)
	{
	  State& invalid_state  = *states[ badStatesIDs[index] ];
	  fout << invalid_state.getCoordinates() << " " << invalid_state << std::endl;
	} // for all invalid states
      
      fout.close();
      
    } // writing badStates
  
  if(isTimeAccurate && isTimeStepTooLarge)
      CFLog(WARN, "The chosen time step is too large as it gives a maximum CFL of " << maxCFL <<".\n");

  if(isTimeAccurate)
      SubSystemStatusStack::getActive()->setMaxDT(SubSystemStatusStack::getActive()->getDT()/maxCFL);
  
  if(isGlobalTimeStep) {
    // AL: I am assuming constant volume here ... needs to be modified with non uniform mesh!!!
    SubSystemStatusStack::getActive()->setMaxDT(CFL/maxUpdateCoeff*volumes[0]);
  }
  
  // computation of the norm of the dU (a.k.a rhs)
  CFreal value = 0.0;
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const CFreal tmp = rhs(iState, getMethodData().getVarID(), nbEqs);
    value += tmp*tmp;
  }
  
  CFreal invalue = value;
  value=0.;
  
  const std::string nsp = getMethodData().getNamespace();
  
  MPI_Allreduce(&invalue,&value,1, Common::MPIStructDef::getMPIType(&value),
		MPI_SUM, PE::GetPE().GetCommunicator(nsp)); // THIS IS A HACK, DO IT PROPERLY!!!
  
  //    CF_DEBUG_OBJ(value);
  value = log10(sqrt(value));
  // CF_DEBUG_OBJ(value);
  getMethodData().setNorm(value);
  //  getMethodData().computeNorm();
  
  //   SafePtr<OutputFormatter> output =
  //     getMethodData().getCollaborator<OutputFormatter>();
  
  //   CFLog(INFO, "Writing output file ... \n");
  //   output->open();
  //   output->write();
  //   output->close();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UpdateSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
   result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSol::setup()
{ 
  const CFuint nbStates = socket_states.getDataHandle().size();
  cf_assert(nbStates > 0);
  m_activeStates.resize(nbStates, true);
  
  CFuint countInactive = 0;
  for (CFuint i =0; i < m_trsWithNullRHS.size(); ++i) {
    const string trsName = m_trsWithNullRHS[i];
    SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(trsName);
    SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
    for (CFuint s = 0; s < trsStates->size(); ++s) {
      m_activeStates[(*trsStates)[s]] = false;
      countInactive++;
    }
  }
  
  if (countInactive > 0) {
    CFLog(INFO, "UpdateSol::setup() => deactivated " << countInactive << " states \n");
  }
}
    
//////////////////////////////////////////////////////////////////////////////
    
  } // namespace ForwardEuler
  
} // namespace COOLFluiD
