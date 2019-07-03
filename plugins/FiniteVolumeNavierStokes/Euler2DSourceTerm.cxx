#include "FiniteVolumeNavierStokes/Euler2DSourceTerm.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Euler2DSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeNavierStokesModule>
euler2DSTFVMCCProvider("Euler2DST");

//////////////////////////////////////////////////////////////////////////////

Euler2DSourceTerm::Euler2DSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _temp(),
  _physicalData()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DSourceTerm::~Euler2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();

  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  cf_assert(_varSet.isNotNull());
  _varSet->getModel()->resizePhysicalData(_physicalData); 
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSourceTerm::computeSource(Framework::GeometricEntity *const element,
				      RealVector& source,
				      RealMatrix& jacobian)
{
  CFLogDebugMin( "Euler2DSourceTerm::computeSource()" << "\n");
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  cf_assert(_varSet.isNotNull());
  State& currState = *socket_states.getDataHandle()[element->getID()];
  
  _varSet->computePhysicalData(currState, _physicalData);
   
  
  const CFreal x = currState.getCoordinates()[XX];
  //std::cout << "x = " << x << endl;
  const CFreal y = currState.getCoordinates()[YY];
  //std::cout << "y = " << y << endl;
  const CFreal z = currState.getCoordinates()[ZZ];
  //std::cout << "z = " << z << endl;
  const CFreal rho = std::sqrt(x*x + y*y);
  //std::cout << "rho = " << rho << endl;
  const CFreal r = std::sqrt(x*x + y*y + z*z);
  //std::cout << "r = " << r << endl;
  const CFreal r2 = r*r;
  //std::cout << "r2 = " << r2 << endl;


  const CFreal density  = _physicalData[EulerTerm::RHO];
  //std::cout << "density = " << density << endl;
  const CFreal u    = _physicalData[EulerTerm::VX];
  //std::cout << "u = " << u << endl;
  const CFreal v    = _physicalData[EulerTerm::VY];
  //std::cout << "v = " << v << endl;
  const CFreal w    = _physicalData[EulerTerm::VZ];
  //std::cout << "w = " << w << endl;
  
  const CFreal GMsun = 1.327474512e20;
  CFreal fx = -GMsun/r2*density*x/r;
  //std::cout << "fx = " << fx << endl;
  CFreal fy = -GMsun/r2*density*y/r;
  //std::cout << "fy = " << fy << endl;
  CFreal fz = -GMsun/r2*density*z/r;
  //std::cout << "fz = " << fz << endl;
  CFreal Vr = x/r*u + y/r*v + z/r*w;
  //std::cout << "Vr = " << Vr << endl;
  CFreal gravEnergy = -density*(GMsun/r2)*Vr;
  //std::cout << "gravEnergy = " << gravEnergy << endl;


  source[0] = 0;
  source[1] = fx;
  source[2] = fy;
  source[3] = fz;
  source[4] = gravEnergy;
  
  // here it is assumed that the symmetry plane is located at y = 0!!!
  
  /// @todo check if this should be with or without sign !!!!!
  const CFreal avRadius = (currState.getCoordinates())[1];
  source *= -1./avRadius*volumes[element->getID()]; 
  CFLog(DEBUG_MAX, "Euler2DSourceTerm::computeSource() => source = " << source << "\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
