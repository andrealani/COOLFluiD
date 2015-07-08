#include "SA3DSourceTerm.hh"
#include "SA3DDES97.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolumeTurb/FiniteVolumeSA.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "NavierStokes/EulerVarSet.hh"

///////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::SA;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SA3DDES97,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeSAModule>
SA3DDES97Provider("SA3DDES97");      
      
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SA3DDES97::SA3DDES97(const std::string& name) :
  SA3DSourceTerm(name),
  socket_length_scale("subgrid"),
  socket_switch_function("switch")
  {
    
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SA3DDES97::~SA3DDES97()
{
}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
SA3DDES97::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = 
	SA3DSourceTerm::providesSockets();
  result.push_back(&socket_length_scale);// a pointer showing the values of the SGS
  result.push_back(&socket_switch_function);// a pointer showing the switch function
  
  return result;
}


//////////////////////////////////////////////////////////////////////////////

void SA3DDES97::setup()
{
  CFAUTOTRACE;
  
  SA3DSourceTerm::setup();
  
  // Get number of states to resize the datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  
  // initialize the subgrid scale
  DataHandle< CFreal> SGS = socket_length_scale.getDataHandle();
  SGS.resize(nbStates);
  SGS = 0.0;

  DataHandle< CFreal> m_switch_function = socket_switch_function.getDataHandle();
  m_switch_function.resize(nbStates);
  m_switch_function = 0.0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal SA3DDES97::getSubLenScale (Framework::GeometricEntity *const element)

{
  ///@todo use the new definition of Spalart for the SGS
  
  /// The SGS is calculated as the square root of the cell surface
  
  //get access to the surfaces of the grid cells
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  //give the ID of the current cell
  const CFuint elemID = element->getID(); // use this
  
  CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  
  const CFreal SubLenScale = pow(volumes[elemID],1./CFreal(dim));
  
  // This is the fundamental constant of the DES models. It may be different according to the application
  ///@see Shur et all 1999
  const CFreal ConstDES = 0.65;
  
  return (ConstDES*SubLenScale);
  
}


///////////////////////////////////////////////////////////////////////////////

CFreal SA3DDES97::getDistance (Framework::GeometricEntity *const element)
{

    const CFreal m_WallDistance = SA3DSourceTerm::getDistance (element);
    
    const CFreal m_SubLenScale = SA3DDES97::getSubLenScale(element);
    
    const CFreal d = min(m_WallDistance , m_SubLenScale);
    
    // Added to print the switch function
    
    //Set the physical data for the cell considered
    State *const currState = element->getState(0); 
    
    DataHandle< CFreal> SGS = socket_length_scale.getDataHandle();
    
    SGS[currState->getLocalID()] = d;
    
    DataHandle< CFreal> m_switch_function = socket_switch_function.getDataHandle();
    
    if (d == m_WallDistance)
    {
      m_switch_function[currState->getLocalID()] = 0.;
    }
    else
    {
      m_switch_function[currState->getLocalID()] = 1.;
    }
    
    
    return d;
      
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

       } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////////////////////////////////

