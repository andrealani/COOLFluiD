#include "SA2DSourceTerm.hh"
#include "SA2DDES97.hh"
#include "SA2DDelDES.hh"
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

MethodStrategyProvider<SA2DDelDES,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeSAModule>
SA2DDelDESProvider("SA2DDelDES");      
      
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SA2DDelDES::SA2DDelDES(const std::string& name) :
  SA2DDES97(name)
  {
    
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SA2DDelDES::~SA2DDelDES()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void SA2DDelDES::setup()
{
  CFAUTOTRACE;
  
  SA2DDES97::setup();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal SA2DDelDES::computeDelayingFunction (Framework::GeometricEntity *const element)
{
  CFreal m_totalViscosity = SA2DSourceTerm::getLamViscosity();
  m_totalViscosity += SA2DSourceTerm::getTurbViscosity();
  
  //calculate the square root of the sum of the gradient velocity
  CFreal m_SquareRootOfVelGrads = std::sqrt(SA2DSourceTerm::compSumOfVelocityGrads());
  
  // to avoid division with 0
  ///@see Shur et all "A hybrid RANS-LES approach with delayed-DES and wall-modelled LES capabilities
  m_SquareRootOfVelGrads = max(m_SquareRootOfVelGrads, 1.e-10);
  
  //take distance to the distance to the wall 
  const CFreal m_WallDistance = SA2DSourceTerm::getDistance(element);
  
  //the Karman Constant
  const CFreal kappa = 0.41;
  
  const CFreal rd = m_totalViscosity/(m_SquareRootOfVelGrads*m_WallDistance*m_WallDistance*kappa*kappa);
  
  const CFreal cubicEightrd = std::pow(8.0*rd , 3);
  
  const CFreal delayingFun = 1. - std::tanh(cubicEightrd);
  
  return delayingFun;
}


///////////////////////////////////////////////////////////////////////////////

CFreal SA2DDelDES::getDistance (Framework::GeometricEntity *const element)
{

    const CFreal m_WallDistance = SA2DSourceTerm::getDistance (element);
    
    const CFreal m_SubLenScale = SA2DDES97::getSubLenScale (element);
    
    const CFreal m_delayingFun = SA2DDelDES::computeDelayingFunction (element);
    
    const CFreal d = m_WallDistance - m_delayingFun*(std::max( 0.0 , (m_WallDistance - m_SubLenScale) ) );
    
    return d;
      
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

       } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////////////////////////////////

