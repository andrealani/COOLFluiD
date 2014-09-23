#include "FiniteVolumeNavierStokes/Euler2DSourceTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
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

  cf_assert(_varSet.isNotNull());
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
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
  
  const CFreal rho  = _physicalData[EulerTerm::RHO];
  const CFreal u    = _physicalData[EulerTerm::VX];
  const CFreal v    = _physicalData[EulerTerm::VY];
  const CFreal rhov = rho*v;
  
  source[0] = rhov;
  source[1] = rhov*u;
  source[2] = rhov*v;
  source[3] = rhov*_physicalData[EulerTerm::H];
  
  // here it is assumed that the symmetry plane is located at y = 0!!!
  
  /// @todo check if this should be with or without sign !!!!!
  const CFreal avRadius = (currState.getCoordinates())[1];
  source *= -1./avRadius*volumes[element->getID()];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
