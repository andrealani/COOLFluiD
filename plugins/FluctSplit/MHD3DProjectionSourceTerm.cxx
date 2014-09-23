#include <numeric>

#include "MHD3DProjectionSourceTerm.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
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

MethodStrategyProvider<MHD3DProjectionSourceTerm,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitMHDModule>
mhd3DProjectionSTProvider("MHD3DProjectionST");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionSourceTerm::MHD3DProjectionSourceTerm(const std::string& name) :
  FluctSplit::ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  m_physicalData()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionSourceTerm::~MHD3DProjectionSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionSourceTerm::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull()); 
  _varSet->getModel()->resizePhysicalData(m_physicalData);
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
			                          RealVector& source,
			                          const FluctSplit::InwardNormalsData& normalsData)
{
  // const vector<State*>& states = *cell->getStates(); 
  
//   CFreal Bdotn = 0.;
//   for (CFuint iState = 0; iState < states.size(); ++iState) {
//     _varSet->computePhysicalData(*states[iState], _physicalData);

//     Bdotn += _physicalData[MHDProjectionTerm::BX] * normalsData.getNodalNormComp(iState,0) +
//       _physicalData[MHDProjectionTerm::BY] * normalsData.getNodalNormComp(iState,1);
//   }

//   Bdotn *= 0.5;
  
//   source = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
