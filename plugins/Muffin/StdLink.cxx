
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/System.hh"
#include "Muffin/CC.hh"
#include "Muffin/BC.hh"
#include "Muffin/StdLink.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdLink,MuffinData,MuffinModule > cStdLinkProvider("StdLink");

//////////////////////////////////////////////////////////////////////////////

void StdLink::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  // get System's to setup, plus CC's and BC's available
  std::vector< SafePtr< System > >& vs = d.m_vcomm_sys;
  std::vector< SafePtr< CC > >& all_cc = d.m_vcomm_cc;
  std::vector< SafePtr< BC > >& all_bc = d.m_vcomm_bc;

  std::vector< SafePtr< CC > >::iterator icc;
  std::vector< SafePtr< BC > >::iterator ibc;
  for (unsigned i=0; i<vs.size(); ++i) {

    // link universal boundary conditions (SystemPLaS is exempt)
    if (!vs[i]->hasTag("Muffin::SystemPLaS")) {

      for (ibc=all_bc.begin(); ibc!=all_bc.end(); ++ibc) {
        if ((*ibc)->hasTag("Muffin::BCFile")               ||
            (*ibc)->hasTag("Muffin::BCFunction")           ||
            (*ibc)->hasTag("Muffin::BCMassFlux"))
          (*ibc)->attachTag("Muffin::System:"+vs[i]->getName());
      }

    }

    // link by tagging appropriate coupling and boundary conditions
    if (vs[i]->hasTag("Muffin::SystemFlow")) {

      for (icc=all_cc.begin(); icc!=all_cc.end(); ++icc) {
        if ((*icc)->hasTag("Muffin::CCBuoyancyFromMITReM") ||
            (*icc)->hasTag("Muffin::CCMomentumFromBubbles"))
          (*icc)->attachTag("Muffin::System:" + vs[i]->getName());
      }
      for (ibc=all_bc.begin(); ibc!=all_bc.end(); ++ibc) {
        if ((*ibc)->hasTag("Muffin::BCElectroosmoticWall") ||
            (*ibc)->hasTag("Muffin::BCFixPressure")        ||
            (*ibc)->hasTag("Muffin::BCFixVelocity")        ||
            (*ibc)->hasTag("Muffin::BCWall"))
          (*ibc)->attachTag("Muffin::System:"+vs[i]->getName());
      }

    }
    else if (vs[i]->hasTag("Muffin::SystemMITReM")) {

      for (ibc=all_bc.begin(); ibc!=all_bc.end(); ++ibc) {
        if ((*ibc)->hasTag("Muffin::BCBulk")      ||
            (*ibc)->hasTag("Muffin::BCElectrode"))
          (*ibc)->attachTag("Muffin::System:"+vs[i]->getName());
      }

    }
    else if (vs[i]->hasTag("Muffin::SystemPLaS")) {


    }
    else if (vs[i]->hasTag("Muffin::SystemTemp")) {

      for (ibc=all_bc.begin(); ibc!=all_bc.end(); ++ibc) {
        if ((*ibc)->hasTag("Muffin::BCWall"))
          (*ibc)->attachTag("Muffin::System:"+vs[i]->getName());
      }

    }
    else if (vs[i]->hasTag("Muffin::SystemTurb")) {

      for (ibc=all_bc.begin(); ibc!=all_bc.end(); ++ibc) {
        if ((*ibc)->hasTag("Muffin::BCWall"))
          (*ibc)->attachTag("Muffin::System:"+vs[i]->getName());
      }

    }
  }

  // display linked coupling conditions
  for (icc=all_cc.begin(); icc!=all_cc.end(); ++icc) {
    std::ostringstream m;
    m << "Linked CC \"" << (*icc)->getName() << "\": ";
    (*icc)->print(m," ");
    d.log(m.str());
  }

  // display linked boundary conditions
  for (ibc=all_bc.begin(); ibc!=all_bc.end(); ++ibc) {
    std::ostringstream m;
    m << "Linked BC \"" << (*ibc)->getName() << "\": ";
    (*ibc)->print(m," ");
    d.log(m.str());
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

