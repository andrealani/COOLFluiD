#include "FluctSplit/SuperInletInterpImpl.hh"
#include "FluctSplit/FluctSplit.hh"

#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Config/PositiveLessThanOne.hh"

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Common/OldLookupTable.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInletInterpImpl, FluctuationSplitData, FluctSplitModule> 
superInletInterpImplProvider("SuperInletInterpImpl");

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterpImpl::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
     ("DiagCoeffFactor","Factor to control the diagonal coefficient.");

  options.addConfigOption< std::string>("InputFileName", "Input data file from which all variables will be extrapolated");

}

//////////////////////////////////////////////////////////////////////////////

SuperInletInterpImpl::SuperInletInterpImpl(const std::string& name) :
  SuperInlet(name),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_updateCoeff("updateCoeff"),
  socket_discardTimeJacob("discardTimeJacob"),
  m_lookupState()
{
   addConfigOptionsTo(this);

   _diagCoeffFactor = 1.0;
   setParameter("DiagCoeffFactor",&_diagCoeffFactor);

  m_infile = "";
  setParameter("InputFileName",&m_infile); 

}

//////////////////////////////////////////////////////////////////////////////

SuperInletInterpImpl::~SuperInletInterpImpl()
{
  for (CFuint i = 0; i < m_lookupState.size(); ++i) {
    deletePtr(m_lookupState[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SuperInletInterpImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = SuperInlet::needsSockets();

  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_discardTimeJacob);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterpImpl::setup() 
{
  SuperInlet::setup();

  std::vector< Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< vector<bool> > discardTimeJacob =
    socket_discardTimeJacob.getDataHandle();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for(CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTRS];
    
    CFLogDebugMax("SuperInletInterpImpl::setup() called for TRS: "
                  << trs->getName() << "\n");
    
    Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
    
    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint localStateID = (*statesIdx)[iState];
      
      /// @TODO AL: here check that the order in which BCs are setup
      /// is the order of their execution !!!
      discardTimeJacob[localStateID].resize(nbEqs);
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
		discardTimeJacob[localStateID][iVar] = true;
      }
    }
  }

  m_variables.resize(PhysicalModelStack::getActive()->getDim());

  const std::string nsp = getMethodData().getNamespace();
  
  if (m_infile != "") {
    if (Common::PE::GetPE().IsParallel()) {
      Common::PE::GetPE().setBarrier(nsp);
      for (CFuint p = 0; p < Common::PE::GetPE().GetProcessorCount(nsp); ++p) {
	if (p == Common::PE::GetPE().GetRank(nsp)) {
	  fillTable();
	}
	Common::PE::GetPE().setBarrier(nsp);
      }
    }
    else {
      fillTable();
    }
  }
  else {
    CFLog(WARN, "WARNING: SuperInletInterp::setup() => filename not specified!\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterpImpl::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "SuperInletInterpImpl::execute() called for TRS: "  << trs->getName() << "\n");

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< std::valarray<Framework::State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  bool isUnsteady(false);
  if(SubSystemStatusStack::getActive()->getDT() > 0.) isUnsteady = true;

  RealVector variables(dim);
  if(isUnsteady) variables.resize(dim+1);

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  // block accumulator 1*1
  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(1, 1, nbEqs));

  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal invcfl = 1.0 / cfl;

  RealVector inletState(nbEqs);

  bool applyBC = true;
  State dimState;
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];

    if (state->isParUpdatable()) {

      // find if we should apply the BC to this node
      if (m_checkCondition) {
        const CFreal applyBCvalue = m_condition.Eval(&state->getCoordinates()[0]);
        applyBC = ((!isUpdated[stateID]) && (applyBCvalue > 0.0));
      }
      else {
        applyBC = (!isUpdated[stateID]);
      }

      // apply the BC to this node
      if (applyBC){
        //Set the values of the variables xyz + time
        for (CFuint i = 0; i < state->getCoordinates().size();++i){
          variables[i] = state->getCoordinates()[i];
        }
        if(isUnsteady) variables[state->getCoordinates().size()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

		/* A */
//         //Evaluate the function
//         m_vFunction.evaluate(variables,*_input);
		/* end - A*/
		
		/* B */
		const RealVector& temp = state->getCoordinates();
		for (CFuint i = 0; i < dim; ++i){
		  m_variables[i] = temp[i];
		}
		
		const CFuint nbEqs  = PhysicalModelStack::getActive()->getNbEq();
		const CFreal yCoord = temp[YY];
		
		for (CFuint i = 0; i < nbEqs; ++i) {
		  // interpolated state value in input variables
		  (*_input)[i] = m_lookupState[i]->get(yCoord);
		}
		/* end - B*/
	  
	  // if some interactive variable IDs are specified, multiply
	  // those variables by the given factor
	  if (_interVarIDs.size() > 0) {
		for (CFuint i = 0; i < _interVarIDs.size(); ++i) {
		  (*_input)[_interVarIDs[i]] *= _interFactor;
		}
	  }		
		
	
        //Set the state value
        if (_inputAdimensionalValues){
          inletState = *_inputToUpdateVar->transform(_input);
        }
        else{
          dimState = *_inputToUpdateVar->transform(_input);
          _varSet->setAdimensionalValues(dimState, inletState);
        }

        // compute coefficients
        const CFreal kcoeff = updateCoeff[stateID]*invcfl;
        const CFreal coeff  = max(kcoeff, _diagCoeffFactor*invcfl);

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          // must include inletState so it will work for unsteady BCs
          // which are updated every iteration
          rhs(stateID, iEq, nbEqs) = coeff * (inletState[iEq] - (*state)[iEq]);

          /// @note AL: the following should be wrong ...
          // (coeff + kcoeff) * (inletState[iEq] - (*state)[iEq]);
        }

        acc->setRowIndex(0, stateID);

        // here we have to know how many and which vertices
        // reference the current boundary node to avoid VERY
        // expensive reallocations!!!!
        const CFuint nbEntries = bStatesNeighbors[stateID].size();
        cf_assert(nbEntries > 0);

        for (CFuint i = 0; i < nbEntries; ++i) {
          const CFuint entryID = bStatesNeighbors[stateID][i]->getLocalID();
          acc->setColIndex(0, entryID);
          acc->setValue(0.0);
          if (entryID == stateID) {
            for (CFuint ib = 0; ib < nbEqs; ++ib) {
              acc->setValue(0,0,ib,ib, coeff);
            }
          }
          jacobMatrix->setValues(*acc);
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
      }
    }
  }

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletInterpImpl::fillTable()
{
  boost::filesystem::path filepath = Environment::DirPaths::getInstance().getWorkingDir() / m_infile;
  
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  std::ifstream& fin = fhandle->open(filepath);
  
  std::string variables;
  // read the first line with the variables names
  getline(fin, variables);
  CFuint nbPoints = 0;
  fin >> nbPoints; std::cout << "nbPoints = " << nbPoints << "\n";
  
  // allocate the look up tables
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_lookupState.resize(nbEqs);
  for (CFuint i = 0; i < nbEqs; ++i) {
    m_lookupState[i] = new Common::LookUpTable<CFreal, CFreal>(nbPoints);
  } 
  
  // nbEqs + "y" AL: hardcoded here
  CFreal ycoord;
  CFreal tmpVar;
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    fin >> ycoord;
    for (CFuint i = 0; i < nbEqs; ++i) {
      fin >> tmpVar;
      m_lookupState[i]->insert(ycoord, tmpVar);
    }
  }
  
  // sort the data 
  for (CFuint i = 0; i < nbEqs; ++i) {
    m_lookupState[i]->sortKeys();
  }
  
  fin.close();
}

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FluctSplit
  
} // namespace COOLFluiD
