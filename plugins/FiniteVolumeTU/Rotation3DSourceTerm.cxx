#include "Rotation3DSourceTerm.hh"
#include "NavierStokes/Euler3DRotationVarSet.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeTU/FiniteVolumeTU.hh"
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

MethodStrategyProvider<Rotation3DSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeTUModule>
rotation3DSTFVMCCProvider("Rotation3DST");

//////////////////////////////////////////////////////////////////////////////

Rotation3DSourceTerm::Rotation3DSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _temp(),
  _physicalData()
{
}

//////////////////////////////////////////////////////////////////////////////

Rotation3DSourceTerm::~Rotation3DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Rotation3DSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler3DRotationVarSet>();
  
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  
  cf_assert(_varSet.isNotNull());
  
  _varSet->getModel()->resizePhysicalData(_physicalData);
}

//////////////////////////////////////////////////////////////////////////////

void Rotation3DSourceTerm::computeSource(Framework::GeometricEntity *const element,
					 RealVector& source,
					 RealMatrix& jacobian)
{
  CFLogDebugMin( "Rotation3DSourceTerm::computeSource()" << "\n");
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  cf_assert(_varSet.isNotNull());

  // this source term is for axisymmetric flows
  const vector<State*>& states = *element->getStates();

  cf_assert(states.size() == 1);
  
  State& currState = *states[0];

  _varSet->computePhysicalData(currState, _physicalData);
  
  const CFreal rhoOmegaV = _varSet->getModel()->getOmega()*
    _physicalData[EulerTerm::RHO]*volumes[element->getID()];
  const CFreal v    = _physicalData[EulerTerm::VY];
  const CFreal w    = _physicalData[EulerTerm::VZ];
  
  source[0] = 0;
  source[1] = 0;
  source[2] = rhoOmegaV*w;
  source[3] = - rhoOmegaV*v;
  source[4] = 0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
