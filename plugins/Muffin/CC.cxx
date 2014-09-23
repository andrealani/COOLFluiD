
#include "Muffin/System.hh"
#include "Muffin/CC.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

CC::CC(const std::string& name) :
    MuffinCom(name),
    s_rhs("rhs"),                     // socket sinks
    s_nodes("nodes"),                 //  ...
    s_states("states"),               //  ...
    s_faceneighcell("faceNeighCell")  //  ...
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  m_applyfrom_str = "";
  setParameter("applyFrom",&m_applyfrom_str);
}

//////////////////////////////////////////////////////////////////////////////

void CC::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;

  options.addConfigOption< std::string >("applyFrom","System name to apply this coupling from (default \"\")");
}

//////////////////////////////////////////////////////////////////////////////

void CC::setup()
{
  CFAUTOTRACE;
  MuffinCom::setup();

  // find System to apply from
  const std::vector< SafePtr< System > >& vs = getMethodData().m_vcomm_sys;
  for (unsigned i=0; i<vs.size() && m_applyfrom.isNull(); ++i)
    if (m_applyfrom_str==vs[i]->getName())
      m_applyfrom = vs[i];
  if (m_applyfrom.isNull())
    err("system name was not found: \"" + m_applyfrom_str + "\"");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

