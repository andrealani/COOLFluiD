
#include "Common/BadValueException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/System.hh"
#include "Muffin/Loop.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

// instantiate providers
MethodCommandProvider< LoopIterate,MuffinData,MuffinModule > cLoopIterate("Iterate");
MethodCommandProvider< LoopTimeStep,MuffinData,MuffinModule > cLoopTimeStep("TimeStep");
MethodCommandProvider< LoopConvergeResidual,MuffinData,MuffinModule > cLoopConvergeResidual("ConvergeResidual");
MethodCommandProvider< LoopConvergeSolution,MuffinData,MuffinModule > cLoopConvergeSolution("ConvergeSolution");

//////////////////////////////////////////////////////////////////////////////

Loop::Loop(const std::string& name) :
    MuffinCom(name),
    m_master(false)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  m_niterations = 1;
  m_clevel      = -8.;
  m_cv.clear();
  m_cei.clear();
  m_tduration = 1.;
  m_tstep     = 1.;
  m_tstart    = -1.;
  setParameter("NbIterations",&m_niterations);
  setParameter("ConvergeL2Level",&m_clevel);
  setParameter("ConvergeVariable",&m_cv);
  setParameter("ConvergeEquation",&m_cei);
  setParameter("TimeDuration",&m_tduration);
  setParameter("TimeStep",&m_tstep);
  setParameter("TimeStart",&m_tstart);

  m_command_names.clear();
  setParameter("CommandNames",&m_command_names);
}

//////////////////////////////////////////////////////////////////////////////

void Loop::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;

  options.addConfigOption< CFuint >("NbIterations","Iterate: number of iterations (default 0)");
  options.addConfigOption< CFreal >("ConvergeL2Level","ConvergeResidual/Solution: logarithm L2 maximum level (default -8.)");
  options.addConfigOption< std::vector< std::string > >("ConvergeVariable","ConvergeSolution: variable(s) to be tested (default < >)");
  options.addConfigOption< std::vector< CFuint > >("ConvergeEquation","ConvergeResidual: equation index (indices) to be tested (default < >)");
  options.addConfigOption< CFreal >("TimeDuration","TimeStep: loop time duration (default 0.)");
  options.addConfigOption< CFreal >("TimeStep","TimeStep: loop time step (default 1.)");
  options.addConfigOption< CFreal >("TimeStart","TimeStep: loop starting time, set by first loop (default negative, don't set)");

  options.addConfigOption< std::vector< std::string > >("CommandNames","Command names to loop on, either Systems or other Loops (default <>)");
}

//////////////////////////////////////////////////////////////////////////////

void Loop::setup()
{
  CFAUTOTRACE;
  MuffinCom::setup();
  MuffinData& d = getMethodData();

  // set MuffinCom pointers to loop on
  if (!m_command_names.size())
    return;
  m_commands.clear();
  const std::vector< SafePtr< System    > >& vs = d.m_vcomm_sys;
  const std::vector< SafePtr< Loop      > >& vl = d.m_vcomm_loops;
  const std::vector< SafePtr< MuffinCom > >& vc = d.m_vcomm_std;
  for (unsigned c=0; c<m_command_names.size(); ++c) {
    bool found = false;

    // look in the System names
    for (unsigned i=0; i<vs.size() && !found; ++i)
      if (m_command_names[c] == vs[i]->getName()) {
        m_commands.push_back(vs[i].d_castTo< MuffinCom >());
        found = true;
      }

    // look in the Loop names
    for (unsigned i=0; i<vl.size() && !found; ++i)
      if (m_command_names[c] == vl[i]->getName()) {
        m_commands.push_back(vl[i].d_castTo< MuffinCom >());
        found = true;
      }

    // look into the Standard commands names
    for (unsigned i=0; i<vc.size() && !found; ++i)
      if (m_command_names[c] == vc[i]->getName()) {
        m_commands.push_back(vc[i]);
        found = true;
      }

    if (!found)
      err("command name was not found: \"" + m_command_names[c] + "\"");
  }

  // get indices of variables to converge on
  m_cvi.clear();
  const std::vector< std::string >& vnames = d.m_varnames;
  for (unsigned i=0; i<m_cv.size() && hasTag("LoopConvergeSolution"); ++i) {
    bool found = false;
    for (unsigned j=0; j<vnames.size() && !found; ++j) {
      if (m_cv[i]==vnames[j]) {
        m_cvi.push_back(j);
        found = true;
      }
    }
    if (!found)
      err("variable was not found: \"" + m_cv[i] + "\"");
  }
  if (hasTag("LoopConvergeSolution") && !m_cvi.size())
    err("variables at \"ConvergeVariable\" not set");

  // get indices of equations to converge on
  if (hasTag("LoopConvergeResidual")) {
    if (!m_cei.size())
      err("equation indices at \"ConvergeEquation\" not set");
    else if (*max_element(m_cei.begin(),m_cei.end())>=vnames.size())
      err("equation maximum index is: " + StringOps::to_str(vnames.size()-1));
  }

  // set global loop properties
  d.m_status->setDTDim(0.);
  if (hasTag("LoopTimeStep") && m_tstart>=0.)
    d.m_status->setCurrentTimeDim(m_tstart);

  // set convergence file, with headers for iteration, time, variable,
  // equation indices and loop name
  m_cvg_file = d.getFilename("convergence.txt");
  const std::string nsp = getMethodData().getNamespace();
  if (!PE::GetPE().GetRank(nsp)) {
    std::ofstream o(m_cvg_file.c_str(),std::ios::trunc);
    o << 'i' << '\t'
      << 't' << '\t'
      << 'n' << '\t';
    for (int e=0; e<Neqns; ++e)
      o << '\"' << d.m_varnames[e] << '\"' << '\t';
    for (int e=0; e<Neqns; ++e)
      o << 'e' << e << '\t';
    o << std::endl;
    o.close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void Loop::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  // set loop properties
  CFuint i = 0;
  CFreal t = 0.;
  if (hasTag("LoopTimeStep") && m_tstart>=0. && m_master)
      d.m_status->setCurrentTimeDim(m_tstart);
  while (!finish(i,t,d.m_logl2_states,d.m_logl2_rhs)) {
    log("step...");

    // display current status information
    if (m_master)
      ver("begin @" + getStatus(d.m_status->getNbIter(),d.m_status->getCurrentTimeDim(),d.m_logl2_states,d.m_logl2_rhs));

    // execute commands sequentially
    std::vector< SafePtr< MuffinCom > >::iterator c;
    if (hasTag("LoopTimeStep")) {
      // make sure this time-step is imposed to allow nested time-stepping
      // (time is updated only from innermost time loop)
      const CFreal tbefore = d.m_status->getCurrentTimeDim();
      for (c=m_commands.begin(); c!=m_commands.end(); ++c) {
        d.m_status->setDTDim(m_tstep);
        ver("call \"" + (*c)->getName() + "\"...");
        (*c)->execute();
        ver("call \"" + (*c)->getName() + "\".");
        d.m_status->setDTDim(m_tstep);
      }
      const CFreal tafter = d.m_status->getCurrentTimeDim();
      if (tbefore==tafter)
        d.m_status->setCurrentTimeDim(tbefore+m_tstep);
    }
    else {
      for (c=m_commands.begin(); c!=m_commands.end(); ++c) {
        ver("call \"" + (*c)->getName() + "\"...");
        (*c)->execute();
        ver("call \"" + (*c)->getName() + "\".");
      }
    }

    // update loop properties (iteration/time)
    ++i;
    t += m_tstep;
    d.m_status->updateNbIter();

    // write convergence and display current status information
    writeConvergence(d.m_status->getNbIter(),d.m_status->getCurrentTimeDim(),d.m_logl2_states,d.m_logl2_rhs);
    if (m_master)
      log("finish @" + getStatus(d.m_status->getNbIter(),d.m_status->getCurrentTimeDim(),d.m_logl2_states,d.m_logl2_rhs));

    log("step.");
  }
}

//////////////////////////////////////////////////////////////////////////////

bool Loop::finish(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs)
{
  const std::string nsp = getMethodData().getNamespace();
  PE::GetPE().setBarrier(nsp);

  // calculate average of convergence norms
  CFreal logl2 = 0.;
  if (hasTag("LoopConvergeResidual")) {
    for (std::vector< CFuint >::const_iterator i=m_cei.begin(); i!=m_cei.end(); ++i)
      logl2 += logl2_rhs[*i]/CFreal(m_cei.size());
  }
  else if (hasTag("LoopConvergeSolution")) {
    for (std::vector< CFuint >::const_iterator i=m_cvi.begin(); i!=m_cvi.end(); ++i)
      logl2 += logl2_states[*i]/CFreal(m_cvi.size());
  }

  // check convergence with appropriate criteria
  int yes = ((hasTag("LoopIterate")?          i>=m_niterations :
             (hasTag("LoopTimeStep")?         t>=m_tduration || MathTools::MathChecks::isZeroWithError(t-m_tduration,1.e-6) :
             (hasTag("LoopConvergeResidual")? (i) && logl2 < m_clevel :
             (hasTag("LoopConvergeSolution")? (i) && logl2 < m_clevel :
                                              true ))))? 1:0);

  // only if all tribal elders approve
  GlobalReduceOperation< GRO_MIN >(&yes,&yes);
  PE::GetPE().setBarrier(nsp);

  return (yes>0);
}

//////////////////////////////////////////////////////////////////////////////

std::string Loop::getStatus(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs)
{
  std::ostringstream msg;
  msg << " i=" << i
      << " t=" << t;
  for (CFuint e=0; e<logl2_states.size(); ++e)
    msg << " L2[\"" << getMethodData().m_varnames[e] << "\"]=" << logl2_states[e];
  for (CFuint e=0; e<logl2_rhs.size(); ++e)
    msg << " L2[e" << e << "]=" << logl2_rhs[e];
  return msg.str();
}

//////////////////////////////////////////////////////////////////////////////

void Loop::writeConvergence(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs)
{
  const std::string nsp = getMethodData().getNamespace();

  PE::GetPE().setBarrier(nsp);
  if (!PE::GetPE().GetRank(nsp)) {
    std::ofstream o(m_cvg_file.c_str(),std::ios::app);
    o << i << '\t'
      << t << '\t'
      << getName() << '\t';
    for (CFuint e=0; e<logl2_states.size(); ++e)
      o << logl2_states[e] << '\t';
    for (CFuint e=0; e<logl2_rhs.size(); ++e)
      o << logl2_rhs[e] << '\t';
    o << std::endl;
    o.close();
  }
  PE::GetPE().setBarrier(nsp);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

