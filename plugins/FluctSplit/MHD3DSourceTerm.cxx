#include "MHD3DSourceTerm.hh"
#include "MHD/MHD3DVarSet.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "InwardNormalsData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitMHD.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHD3DSourceTerm,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitMHDModule>
mhd3DSTProvider("MHD3DST");

//////////////////////////////////////////////////////////////////////////////

MHD3DSourceTerm::MHD3DSourceTerm(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _physicalData()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DSourceTerm::~MHD3DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DSourceTerm::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();
  cf_assert(_varSet.isNotNull());
  _varSet->getModel()->resizePhysicalData(_physicalData);
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
			                          RealVector& source,
			                          const FluctSplit::InwardNormalsData& normalsData)
{
  const vector<State*>& states = *cell->getStates(); 
  
  const RealVector& linearData = _varSet->getModel()->getPhysicalData();

  const CFreal ubar = linearData[MHDTerm::VX];
  const CFreal vbar = linearData[MHDTerm::VY];
  const CFreal wbar = linearData[MHDTerm::VZ];
  const CFreal Bxbar = linearData[MHDTerm::BX];
  const CFreal Bybar = linearData[MHDTerm::BY];
  const CFreal Bzbar = linearData[MHDTerm::BZ];

  CFreal Bdotn = 0.;
  for (CFuint iState = 0; iState < states.size(); ++iState) {
    _varSet->computePhysicalData(*states[iState], _physicalData);

    Bdotn += _physicalData[MHDTerm::BX] * normalsData.getNodalNormComp(iState,0) +
      _physicalData[MHDTerm::BY] * normalsData.getNodalNormComp(iState,1) +
      _physicalData[MHDTerm::BZ] * normalsData.getNodalNormComp(iState,2);
  }
  
  //Bdotn *= 0.5;
  Bdotn *= (-1./3.);
  
  const CFreal VdotB = ubar*Bxbar + vbar*Bybar + wbar*Bzbar;

  source[1] =  Bxbar*Bdotn;
  source[2] =  Bybar*Bdotn;
  source[3] =  Bzbar*Bdotn;
  source[4] =  ubar*Bdotn;
  source[5] =  vbar*Bdotn;
  source[6] =  wbar*Bdotn;
  source[7] =  VdotB*Bdotn;
  
  //FOR_VERTEX
  //{
  //node[local_node[v]].divB += over_VE*Bdotn;
  //}
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
