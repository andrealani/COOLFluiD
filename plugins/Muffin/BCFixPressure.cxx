
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/System.hh"
#include "Muffin/SystemTurb.hh"
#include "Muffin/BCFixPressure.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCFixPressure,MuffinData,MuffinModule > cBCFixPressureProvider("FixPressure");

//////////////////////////////////////////////////////////////////////////////

BCFixPressure::BCFixPressure(const std::string& name) :
    BC(name)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCFixPressure");
  addConfigOptionsTo(this);

  // set parameters
  m_value = 0.;
  setParameter("Value",&m_value);
}

//////////////////////////////////////////////////////////////////////////////

BCFixPressure::~BCFixPressure()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCFixPressure::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< CFreal >("Value","Value (default 0.)");
}

//////////////////////////////////////////////////////////////////////////////

void BCFixPressure::setup()
{
  CFAUTOTRACE;
  BC::setup();

  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  // for all the boundary regions, cycle the nodes
  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<trs.size(); ++i) {
    const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();
    for (CFuint inb=0; inb<nodes.size(); ++inb)
      setInitialState(nodes[inb],VPRESSURE,m_value);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFixPressure::applyOnSystemFlow(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  const int ik = (turmod? getMethodData().getVariableIndex(VTURBK) : -1);

  // cycle the nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {
    setDirichletCondition( *s,*n,0,
      (*h_states[*n])[s->iv] - (ik<0? 0. : 0.66667*(*h_states[*n])[ik]) );

  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

