#include "SA3DSourceTerm.hh"
#include "SA3DDES97.hh"
#include "SA3DDelDES.hh"
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

MethodStrategyProvider<SA3DDelDES,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeSAModule>
SA3DDelDESProvider("SA3DDelDES");      
      
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SA3DDelDES::SA3DDelDES(const std::string& name) :
  SA3DDES97(name)
  {
    
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SA3DDelDES::~SA3DDelDES()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void SA3DDelDES::setup()
{
  CFAUTOTRACE;
  
  SA3DDES97::setup();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal SA3DDelDES::computeDelayingFunction (Framework::GeometricEntity *const element)
{
  CFreal m_totalViscosity = SA3DSourceTerm::getLamViscosity();
  m_totalViscosity += SA3DSourceTerm::getTurbViscosity();
  
  //calculate the square root of the sum of the gradient velocity
  CFreal m_SquareRootOfVelGrads = std::sqrt(SA3DSourceTerm::compSumOfVelocityGrads());
  
  // to avoid division with 0
  ///@see Shur et all "A hybrid RANS-LES approach with delayed-DES and wall-modelled LES capabilities
  m_SquareRootOfVelGrads = max(m_SquareRootOfVelGrads, 1.e-10);
  
  //take distance to the distance to the wall 
  const CFreal m_WallDistance = SA3DSourceTerm::getDistance(element);
  
  //the Karman Constant
  const CFreal kappa = 0.41;
  
  const CFreal rd = m_totalViscosity/(m_SquareRootOfVelGrads*m_WallDistance*m_WallDistance*kappa*kappa);
  
  const CFreal cubicEightrd = std::pow(8.0*rd , 3);
  
  const CFreal delayingFun = 1. - std::tanh(cubicEightrd);
  
  return delayingFun;
}


///////////////////////////////////////////////////////////////////////////////

CFreal SA3DDelDES::getDistance (Framework::GeometricEntity *const element)
{

    const CFreal m_WallDistance = SA3DSourceTerm::getDistance (element);
    
    const CFreal m_SubLenScale = SA3DDES97::getSubLenScale (element);
    
    const CFreal m_delayingFun = SA3DDelDES::computeDelayingFunction (element);
    
    const CFreal d = m_WallDistance - m_delayingFun*(std::max( 0.0 , (m_WallDistance - m_SubLenScale) ) );
    
    // Added to print the switch function
    
    //Set the physical data for the cell considered
    State *const currState = element->getState(0); 
    
    DataHandle< CFreal> SGS = socket_length_scale.getDataHandle();
    
    SGS[currState->getLocalID()] = d;
    
    DataHandle< CFreal> m_switch_function = socket_switch_function.getDataHandle();
    
      m_switch_function[currState->getLocalID()] = m_delayingFun;
    
    
    return d;
      
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

       } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////////////////////////////////

