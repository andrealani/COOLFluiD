
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/LoopBalanceCurrent.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

// instantiate provider
MethodCommandProvider< LoopBalanceCurrent,MuffinData,MuffinModule > cLoopBalanceCurrent("BalanceCurrent");

//////////////////////////////////////////////////////////////////////////////

LoopBalanceCurrent::LoopBalanceCurrent(const std::string& name) :
    Loop(name)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("LoopConvergeSolution"); // make this work as a LoopConvergeSolution
  attachTag("LoopBalanceCurrent");
  addConfigOptionsTo(this);

  m_current_l2_level = -6.;
  m_sysmitrem_str.clear();
  setParameter("BalanceCurrentL2Level",&m_current_l2_level);
  setParameter("BalanceCurrentSystem",&m_sysmitrem_str);
}

//////////////////////////////////////////////////////////////////////////////

void LoopBalanceCurrent::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< CFreal >("BalanceCurrentL2Level","BalanceCurrent: logarithm L2 of total current maximum level (default -6.)");
  options.addConfigOption< std::string >("BalanceCurrentSystem","MITReM system to access boundary currents from, name (default \"\")");
}

//////////////////////////////////////////////////////////////////////////////

void LoopBalanceCurrent::setup()
{
  CFAUTOTRACE;
  Loop::setup();

  // set pointer to system MITReM, looking in the System names
  // (an empty name will match the first found system)
  const std::vector< SafePtr< System > >& vs = getMethodData().m_vcomm_sys;
  bool found = false;
  for (unsigned i=0; i<vs.size() && !found; ++i) {
    if (vs[i]->hasTag("Muffin::SystemMITReM") &&
       (m_sysmitrem_str.length()? vs[i]->getName()==m_sysmitrem_str : true)) {

      try {
        m_sysmitrem = vs[i].d_castTo< SystemMITReM >();
        found = true;
      }
      catch (FailedCastException& e) {
        log("couldn't cast \"" + vs[i]->getName() + "\" to SystemMITReM");
        throw;
      }

    }
  }
  if (!found) {
    if (m_sysmitrem_str.length())
      err("no MITReM system with name \"" + m_sysmitrem_str + "\" was found");
    else
      err("no MITReM system available");
  }
}

//////////////////////////////////////////////////////////////////////////////

bool LoopBalanceCurrent::finish(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs)
{
  // update current at each TRS's and over all ranks
  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");
  const CFreal eps = MathTools::MathConsts::CFrealEps();
  std::vector< CFreal > vi = m_sysmitrem->m_vi;
  m_current_bal = 0.;
  m_current_max = 0.;
  for (CFuint j=0; j<vi.size(); ++j) {
    GlobalReduceOperation< GRO_SUM >(&vi[j],&vi[j]);
    m_current_bal += vi[j];
    m_current_max = std::max(m_current_max,fabs(vi[j]));
    log("Current at TRS \"" + btrs[j]->getName() + "\" [A]: " + StringOps::to_str(vi[j]));
  }
  m_current_bal = fabs(m_current_bal)+eps;
  m_current_max = std::max(m_current_max,eps);

  const CFreal logl2 = log10(m_current_bal/m_current_max);
  log("Current L2 level: " + StringOps::to_str(logl2));
  return Loop::finish(i,t,logl2_states,logl2_rhs) &&
         (logl2<m_current_l2_level);
}

//////////////////////////////////////////////////////////////////////////////

std::string LoopBalanceCurrent::getStatus(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs)
{
  std::ostringstream msg;
  const CFreal l2 = log10(m_current_bal/m_current_max);
  msg << Loop::getStatus(i,t,logl2_states,logl2_rhs) << " L2[current]=" << l2;
  return msg.str();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

