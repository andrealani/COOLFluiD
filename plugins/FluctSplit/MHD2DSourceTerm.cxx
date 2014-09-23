#include "MHD2DSourceTerm.hh"
#include "MHD/MHD2DVarSet.hh"
#include "Common/CFLog.hh"
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

MethodStrategyProvider<MHD2DSourceTerm,
		       FluctuationSplitData,
		       ComputeSourceTerm<FluctuationSplitData>,
		       FluctSplitMHDModule>
mhd2DSTProvider("MHD2DST");

//////////////////////////////////////////////////////////////////////////////

MHD2DSourceTerm::MHD2DSourceTerm(const std::string& name) :
  FluctSplit::ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  m_physicalData()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DSourceTerm::~MHD2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DSourceTerm::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD2DVarSet>();
  cf_assert(_varSet.isNotNull());
  _varSet->getModel()->resizePhysicalData(m_physicalData); 
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
			                          RealVector& source,
			                          const FluctSplit::InwardNormalsData& normalsData)
{
  const vector<State*>& states = *cell->getStates(); 
    
  const RealVector& linearData =
    _varSet->getModel()->getPhysicalData();

  const CFreal ubar = linearData[MHDTerm::VX];
  const CFreal vbar = linearData[MHDTerm::VY];
  const CFreal wbar = linearData[MHDTerm::VZ];
  const CFreal Bxbar = linearData[MHDTerm::BX];
  const CFreal Bybar = linearData[MHDTerm::BY];
  const CFreal Bzbar = linearData[MHDTerm::BZ];

  CFreal Bdotn = 0.;
  for (CFuint iState = 0; iState < states.size(); ++iState) {
    _varSet->computePhysicalData(*states[iState], m_physicalData);

    Bdotn += m_physicalData[MHDTerm::BX] * normalsData.getNodalNormComp(iState,0) +
      m_physicalData[MHDTerm::BY] * normalsData.getNodalNormComp(iState,1);
  }
  
  Bdotn *= (-0.5);
  
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
