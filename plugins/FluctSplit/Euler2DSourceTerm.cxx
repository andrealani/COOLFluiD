#include <numeric>

#include "Euler2DSourceTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "InwardNormalsData.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Euler2DSourceTerm,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitNavierStokesModule>
euler2DSTProvider("Euler2DAxiST");

//////////////////////////////////////////////////////////////////////////////

Euler2DSourceTerm::Euler2DSourceTerm(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DSourceTerm::~Euler2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSourceTerm::setup()
{
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
  const vector<State*>& states = *cell->getStates(); 
    
  // this source term is for axisymmetric flows
  const RealVector& linearData =
    _varSet->getModel()->getPhysicalData();

  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal rhov = rho*v;
  source[0] = rhov;
  source[1] = rhov*u;
  source[2] = rhov*v;
  source[3] = rhov*linearData[EulerTerm::H];
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  CFreal avRadius = 0.;
  const CFuint nbStatesInCell = states.size();
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    avRadius += (states[iState]->getCoordinates())[1];
  }
  avRadius /= nbStatesInCell;
  
  source *= (-volumes[cell->getID()]/avRadius);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
