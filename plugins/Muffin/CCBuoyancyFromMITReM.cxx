
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/SystemMITReM.hh"
#include "Muffin/SystemFlow.hh"
#include "Muffin/CCBuoyancyFromMITReM.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< CCBuoyancyFromMITReM,MuffinData,MuffinModule > cCCBuoyancyFromMITReM("BuoyancyFromMITReM");

//////////////////////////////////////////////////////////////////////////////

CCBuoyancyFromMITReM::CCBuoyancyFromMITReM(const std::string& name) :
    CC(name)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::CCBuoyancyFromMITReM");
  addConfigOptionsTo(this);

//  m_system_from = "";
//  m_system_to = "";
//  setParameter("SystemFrom",&m_system_from);
//  setParameter("SystemTo",&m_system_to);
}

//////////////////////////////////////////////////////////////////////////////

void CCBuoyancyFromMITReM::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;

//  options.addConfigOption< std::string >("SystemFrom","System to couple from, default \"\"");
//  options.addConfigOption< std::string >("SystemTo","System to couple to, default \"\"");
}

//////////////////////////////////////////////////////////////////////////////

void CCBuoyancyFromMITReM::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  CC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CCBuoyancyFromMITReM::setup()
{
  CFAUTOTRACE;
  CC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void CCBuoyancyFromMITReM::unsetup()
{
  CFAUTOTRACE;
  CC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CCBuoyancyFromMITReM::apply(const Common::SafePtr< System > applyto)
{
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

