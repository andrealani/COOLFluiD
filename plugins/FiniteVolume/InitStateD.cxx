#include "FiniteVolume/FiniteVolume.hh"

#include "InitStateD.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateD, CellCenterFVMData, FiniteVolumeModule>
initStateDProvider("InitStateD");

//////////////////////////////////////////////////////////////////////////////

InitStateD::InitStateD(const std::string& name) :
  InitState(name),
  socket_wallDistance("wallDistance")
{
}

//////////////////////////////////////////////////////////////////////////////

InitStateD::~InitStateD()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitStateD::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitStateD::executeOnTrs() called for TRS: " << trs->getName() << "\n");
  
  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitStateD not applied to InnerFaces!!!");
  }
  
  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector coordD(dim+1);
  State dimState;
  
  std::vector<CFuint>::iterator itd;
  for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
    State* const currState = states[(*itd)];
    coordD.slice(0,dim) = currState->getCoordinates(); 
    cf_assert(currState->getLocalID() < wallDistance.size());
    coordD[dim] = wallDistance[currState->getLocalID()];
    _vFunction.evaluate(coordD, *_input);
    if (_inputAdimensionalValues) { cout << "e" << endl;
      *currState = *_inputToUpdateVar->transform(_input);
    }
    else {
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > InitStateD::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = InitState::needsSockets();
  result.push_back(&socket_wallDistance);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume
    
  } // namespace Numerics
  
} // namespace COOLFluiD
