#include "Common/PE.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/CoupledDirichletBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledDirichletBC, FiniteElementMethodData, FiniteElementModule> CoupledDirichletBCProvider("CoupledDirichletBC");

//////////////////////////////////////////////////////////////////////////////

void CoupledDirichletBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< bool >("Implicit","Apply the BC implicitly?");
   options.addConfigOption< std::string >( "Symmetry",
     "Keep matrix symmetry by AdjustColumn or ScaleDiagonal methods or not (default)." );
   options.addConfigOption< std::vector< CFuint > >( "ApplyEqs",
     "Apply the BC only to specified equations zero-based indexed (default all equations)." );
   options.addConfigOption< CFreal >( "ScaleDiagonal",
     "ScaleDiagonal symmetry method coefficient (default 1.e20)." );
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
   options.addConfigOption< bool >("UseDeltaStates","If true, use (states-pastStates) instead of (states)");
   options.addConfigOption< CFuint >("SubIterations","Number of subiterations done for each otherSubSysStep");

   options.addConfigOption< bool >("AlternateBC","Alternate between different BCs");
//    options.addConfigOption< CFuint >("AlternateRate","Alternate Rate between different BCs");
   options.addConfigOption< bool >("AlternateStart","Start with this BC when alternating between different BCs");
}

//////////////////////////////////////////////////////////////////////////////

CoupledDirichletBC::CoupledDirichletBC(const std::string& name) :
  FiniteElementMethodCom(name),
  _sockets(),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_states("states"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_appliedStrongBC("appliedStrongBC")
{
   addConfigOptionsTo(this);

  m_applyEqs = std::vector< CFuint >();
  m_symmetryStr = "No";
  m_scale = 1.e20;

  _interfaceName = "";
   setParameter("Interface",&_interfaceName);

  _isImplicit = false;
   setParameter("Implicit",&_isImplicit);

  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);

  _useDeltaStates = false;
  setParameter("UseDeltaStates",&_useDeltaStates);

  _nbSubIterations = 1;
  setParameter("SubIterations",&_nbSubIterations);

  _alternateBC = false;
  setParameter("AlternateBC",&_alternateBC);

//   _alternateRate = 2;
//   setParameter("AlternateRate",&_alternateRate);

  _alternateStart = true;
  setParameter("AlternateStart",&_alternateStart);


}

//////////////////////////////////////////////////////////////////////////////

CoupledDirichletBC::~CoupledDirichletBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledDirichletBC::setup()
{
  // first call parent method
  FiniteElementMethodCom::setup();

  _currentSubIteration = 0;
  _previousIteration = 0;

  _currentAlternateRun = _alternateStart;
  if(!_alternateBC) _currentAlternateRun = true;

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
        throw BadValueException (FromHere(),
          "CoupledDirichletBC ApplyEqs refers to an equation that doesn't exist." );
  }

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  DataHandle< vector<bool> > appliedStrongBC =
    socket_appliedStrongBC.getDataHandle();

  std::vector< Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();

  for(CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTRS];

    CFLogDebugMax("CoupledDirichletBC::setup() called for TRS: "
		  << trs->getName() << "\n");

    Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint localStateID = (*statesIdx)[iState];

      appliedStrongBC[localStateID].resize(nbEqs);
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	      appliedStrongBC[localStateID][iVar] =
            (m_applyEqs[iVar]) ? true : false;
      }
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

void CoupledDirichletBC::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledDirichletBC::executeOnTrs()
{
  CFAUTOTRACE;
///@todo HERE subiteration is meant when you have subsys1
///doing 3 steps while subsys2 does 1 step for example!!!
///Change this name!!!

  const CFuint currentIteration = SubSystemStatusStack::getActive()->getNbIter();
//    std::cout << "Iteration: " << currentIteration << std::endl;
//    std::cout << "Previous Iteration: " << _previousIteration << std::endl;
//    std::cout << "SubIteration: " << _currentSubIteration << "/" << _nbSubIterations << std::endl;

  if((currentIteration - _previousIteration) > 0)
  {
    cf_assert((currentIteration - _previousIteration) == 1);
    _currentSubIteration++;
    if(_currentSubIteration > _nbSubIterations) _currentSubIteration -= _nbSubIterations;
  }

_currentSubIteration = 1;
//    std::cout << "SubIteration: " << _currentSubIteration << "/" << _nbSubIterations << std::endl;
  if(_currentAlternateRun)
  {

   std::cout << "Running DirichletBC" << std::endl;

    getMethodData().setDirichletBCApplied(true);

    SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
    CFLogDebugMin( "CoupledDirichletBC::executeOnTrs() called for TRS: " << trs->getName() << "\n");

    DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
    DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
    DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
    DataHandle<std::valarray<Framework::State*> > bStatesNeighbors =
      socket_bStatesNeighbors.getDataHandle();

    SafePtr<LinearSystemSolver> lss =
      getMethodData().getLinearSystemSolver()[0];
    const LSSIdxMapping& idxMapping = lss->getLocalToGlobalMapping();

    SafePtr<LSSMatrix> sysMat = lss->getMatrix();

    // this should be an intermediate lightweight assembly !!!
    // it is needed because here you SET values while elsewhere
    // you ADD values
    sysMat->flushAssembly();

    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint dim   = PhysicalModelStack::getActive()->getDim();
    const CFuint nbApplyEqs = m_applyEqs.size();
    const bool isAdjust = (m_symmetryStr=="AdjustColumn");
    const bool isScale  = (m_symmetryStr=="ScaleDiagonal");
    if (!isScale)
      m_scale = 1.;
    CFreal implicit = (_isImplicit? 1.:0.);

    // block accumulator 1*1
    //auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(1, 1, nbEqs));

    const std::string trsName = trs->getName();
    SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
    std::vector<CFuint>::iterator itd;

    const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
    const std::string nameSpace = getMethodData().getNamespace();

    const std::string baseName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_States_";
    const std::string socketDataName = baseName + "DATA";
    const std::string socketAcceptName = baseName + "ISACCEPTED";

    RealVector dirichletState(nbEqs);
    RealVector dirichletRhs(nbApplyEqs);
    RealVector variables(dim + 1 + nbEqs);

    CFuint idxData(0);
    CFuint idxAccepted(0);
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State *const currState = states[*itd];

      if (currState->isParUpdatable()) {

        const CFuint localStateID = currState->getLocalID();
        const RealVector& temp = currState->getCoordinates() ;
        for (CFuint iCoord = 0; iCoord < temp.size();++iCoord){
          variables[iCoord] = temp[iCoord];
        }
        variables[dim] = SubSystemStatusStack::getActive()->getCurrentTimeDim();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
          variables[dim + 1 + iEq] = (*currState)[iEq];
        }

        if (!isUpdated[localStateID]) {

          DataHandle<RealVector>
            interfaceData = _sockets.getSocketSink<RealVector>(socketDataName)->getDataHandle();
          DataHandle<RealVector>
            interfacePastData = _sockets.getSocketSink<RealVector>(socketDataName + "_PAST")->getDataHandle();
          DataHandle<CFreal>
            isAccepted = _sockets.getSocketSink<CFreal>(socketAcceptName)->getDataHandle();

          if(isAccepted[idxAccepted] >= 0.)
          {
            dirichletState = interfaceData[idxData];
// CFout << "The state: " << currState->getCoordinates() << " has received value: "<< dirichletState << "\n";
            if((_useDeltaStates) &&(SubSystemStatusStack::getActive()->getNbIter() == 1) && SubSystemStatusStack::getActive()->isSubIterationFirstStep())
            {
              dirichletState = 0.;
              interfacePastData[idxData] = interfaceData[idxData];
            }
            else{
              if(_useDeltaStates){
// CFout << "The past data is: " << interfacePastData[idxData] << "\n";
                dirichletState -= interfacePastData[idxData];
                dirichletState /= _nbSubIterations;
              }
            }
// CFout << "Final Dirichlet state: " << dirichletState<< "\n";
          }
          else{
	   CFLogDebugMin("This state has not been accepted!!! : " << currState->getCoordinates() << "\n");
	   CFLogDebugMin("setting it to the default value!!!\n");
            _vFunction.evaluate(variables,dirichletState);
          }

// CFout << "The state: " << currState->getCoordinates() << " is set to: "<< dirichletState << "\n";
          // global position of node, which must have at least one neighbour
          const CFuint nGlobalID = idxMapping.getColID(localStateID)*nbEqs;
          const CFuint nbNeigh = bStatesNeighbors[localStateID].size();
          cf_assert(nbNeigh>0);

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
                bStatesNeighbors[localStateID][j]->getLocalID() )*nbEqs;

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
              
              CFint * iM = new CFint[M];
              CFint * iN = new CFint[N];
              
              for (CFuint i=0; i<nbNeigh; ++i) {
                const CFuint iGlobalID = idxMapping.getColID(
                  bStatesNeighbors[localStateID][i]->getLocalID() )*nbEqs;
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
                  bStatesNeighbors[localStateID][i]->getLocalID();
                for (CFuint j=0; j<nbApplyEqs; ++j) {
                  const CFuint jEq = m_applyEqs[j];
                  for (CFuint k=0; k<nbApplyEqs; ++k)
                    rhs(neighStateID,jEq,nbEqs) -=
                      sysMatM(i*nbApplyEqs+j,k) * dirichletRhs[k];
                }
              }
              sysMatM = 0.;
              sysMat->setValues( M,iM, N,iN, &sysMatM[0] );

              deletePtrArray(iM);
              deletePtrArray(iN);

            } // isAdjust?
          } // !isScale?


          // set rhs and system matrix diagonal terms (scaled)
          for (CFuint j=0; j<nbApplyEqs; ++j) {
            const CFuint jEq=m_applyEqs[j];
            sysMat->setValue(nGlobalID+jEq,nGlobalID+jEq, m_scale);
            rhs(localStateID,jEq,nbEqs) = dirichletRhs[j] * m_scale;
          }

          // flagging is important!
          isUpdated[localStateID] = true;



/*          // coefficient i,i in global linear system
          const CFreal coeff = 1.;

          // rhs is set to
          if(_isImplicit) {
            for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
              rhs(localStateID, iEq, nbEqs) = dirichletState[iEq] - (*currState)[iEq];
            }
          }
          else {
            for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
              rhs(localStateID, iEq, nbEqs) = dirichletState[iEq];
            }
          }

          acc->setRowIndex(0, localStateID);

          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);

          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID = bStatesNeighbors[localStateID][i]->getLocalID();
            acc->setColIndex(0, entryID);
            acc->setValue(0.0);
            if (entryID == localStateID) {
              for (CFuint ib = 0; ib < nbEqs; ++ib) {
                acc->setValue(0,0,ib,ib, coeff);
              }
            }
            jacobMatrix->setValues(*acc);
          }
          isUpdated[localStateID] = true; // flagging is important!!!!!*/
        }
      }
      //We have to increase the counter even if the state was already
      //updated, otherwise, we could assign values to the wrong nodes!!
      DataHandle<CFreal>
        isAccepted = _sockets.getSocketSink<CFreal>(socketAcceptName)->getDataHandle();

      if(isAccepted[idxAccepted] >= 0.){
        idxData++;
      }

      idxAccepted++;
    }

    // this should be an intermediate lightweight assembly !!!
    // it is needed because here you SET values while elsewhere
    // you ADD values
    sysMat->flushAssembly();

  }

  //update currentAlternateRun
///@todo not ok for global subiteration
  _previousIteration = currentIteration;
  if((currentIteration > 0)&&(_alternateBC)) _currentAlternateRun = !_currentAlternateRun;

}

//////////////////////////////////////////////////////////////////////////////

void CoupledDirichletBC::configure ( Config::ConfigArgs& args )
{
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

  const std::string nameSpace = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];

    const std::string baseName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_States_";
    const std::string socketDataBaseName = baseName + "DATA";
    const std::string socketAcceptBaseName = baseName + "ISACCEPTED";

    _sockets.createSocketSink<CFreal>(socketAcceptBaseName);
    _sockets.createSocketSink<RealVector>(socketDataBaseName);
    _sockets.createSocketSink<RealVector>(socketDataBaseName + "_PAST");

  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledDirichletBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();

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
