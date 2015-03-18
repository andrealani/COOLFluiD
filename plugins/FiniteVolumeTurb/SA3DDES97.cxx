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
  SA3DSourceTerm(name)
  {
    
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////

SA3DDES97::~SA3DDES97()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void SA3DDES97::setup()
{
  CFAUTOTRACE;
  
  SA3DSourceTerm::setup();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal SA3DDES97::getSubLenScale (Framework::GeometricEntity *const element)

{
  ///@todo use the new definition of Spalart for the SGS
  
  /// The SGS is calculated as the square root of the cell surface
  
  //get access to the surfaces of the grid cells
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  //give the ID of the current cell
  const CFuint elemID = element->getID();
  
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
    
    return d;
      
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

       } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////////////////////////////////

