
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/SystemPLaS.hh"
#include "Muffin/SystemFlow.hh"
#include "Muffin/CCMomentumFromBubbles.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< CCMomentumFromBubbles,MuffinData,MuffinModule > cCCMomentumFromBubbles("MomentumFromBubbles");

//////////////////////////////////////////////////////////////////////////////

CCMomentumFromBubbles::CCMomentumFromBubbles(const std::string& name) :
    CC(name),
    s_mn_volume("NodalVolume")  // socket sinks
{
  CFAUTOTRACE;

  attachTag("Muffin::CCMomentumFromBubbles");
  addConfigOptionsTo(this);  // calls defineConfigOptions (necessary?)
}

//////////////////////////////////////////////////////////////////////////////

void CCMomentumFromBubbles::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;

//  options.addConfigOption< std::string >("SystemFrom","System to couple from, default \"\"");
//  options.addConfigOption< std::string >("SystemTo","System to couple to, default \"\"");
}

//////////////////////////////////////////////////////////////////////////////

void CCMomentumFromBubbles::apply(const SafePtr< System > applyto)
{
  // get the systems in their appropriate type
  SafePtr< SystemPLaS > sys_plas;
  SafePtr< SystemFlow > sys_flow;
  try {
    sys_plas = m_applyfrom.d_castTo< SystemPLaS >();
    sys_flow = applyto.d_castTo< SystemFlow >();
  }
  catch (Common::FailedCastException& e) {
    log("System \"" + m_applyfrom->getName() + "\" must be of type PLaS!");
    log("System \"" + applyto->getName() + "\" must be of type Flow!");
    throw;
  }

  // get DataHandles and PLaSTracking phase data
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< CFreal > h_mn_volume = s_mn_volume.getDataHandle();
  const PLAS_PHASE_DATA* pd = sys_plas->getPhaseData();

  // add bubble momentum contribution (first comes pressure, then velocities)
  const std::vector< SafePtr< TopologicalRegionSet > > vtrs = sys_flow->getTrsList();
  for (std::vector< SafePtr< TopologicalRegionSet > >::const_iterator
       t=vtrs.begin(); t!=vtrs.end(); ++t) {
    const std::vector< CFuint > vnodes = *(*t)->getNodesInTrs();
    for (std::vector< CFuint >::const_iterator
         n=vnodes.begin(); n!=vnodes.end(); ++n) {
      for (int e=0; e<sys_flow->Nsys; ++e)
        h_rhs(*n,sys_flow->iv+e,Neqns) += pd[*n].dispForce[e]*h_mn_volume[*n];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

