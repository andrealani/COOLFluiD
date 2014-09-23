#include <numeric>

#include "LinearizedEuler3DSourceTermValiant.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplitLinEuler.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler3DSourceTermValiant,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler3DSTProviderValiant("LinEuler3DSourceValiant");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler3DSourceTermValiant::LinearizedEuler3DSourceTermValiant(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp()
{
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler3DSourceTermValiant::~LinearizedEuler3DSourceTermValiant()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceTermValiant::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler3DVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceTermValiant::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
//IMPORTANT!!!!!
// You should not multiply here by delta t, this is done later for the whole fluctuation!!!!!

  const RealVector centroid = cell->computeCentroid();
  DistributionData& ddata = getMethodData().getDistributionData();
  const CFreal time = ddata.time;
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  vector<State*>& states = *ddata.states;
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal pi=3.14159265;
  const CFreal alpha = 5.522;
  const CFreal f = 138.8;
  const CFreal om = 2.*pi*f;
  const CFreal deltar0 = 0.15;
  const CFreal xs = 0.;
  const CFreal ys = 0.;
  const CFreal zs = 0.;

  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;
  source[4] = 0.0;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)  {

    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    const CFreal z = (states[iState]->getCoordinates())[ZZ];
    const CFreal posseidon = cos(0.5*pi*sqrt((x-xs)*(x-xs)+(y-ys)*(y-ys)+(z-zs)*(z-zs))/deltar0);
    const CFreal psi = sin(2.*pi*f*time);
    const CFreal source_vol = 4./3.*pi*deltar0*deltar0*deltar0;
    CFreal function=alpha/source_vol*posseidon*posseidon*psi;

    source[0] += function;
    source[1] += 0.;
    source[2] += 0.;
    source[3] += 0.;
    source[4] += c*c*function;
  }

  const CFuint dimSource = source.size();

  for (CFuint iSource = 0; iSource < dimSource; ++iSource) {
    source[iSource] /= nbStatesInCell;
    source[iSource] *= volumes[cellID]*dt;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
