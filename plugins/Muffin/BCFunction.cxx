
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/System.hh"
#include "Muffin/BCFunction.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCFunction,MuffinData,MuffinModule > cBCFunctionProvider("Function");

//////////////////////////////////////////////////////////////////////////////

BCFunction::BCFunction(const std::string& name) :
    BC(name)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCFunction");
  addConfigOptionsTo(this);

  // set parameters
  m_function_str.clear();
  m_applyvars_str.clear();
  setParameter("Def",&m_function_str);
  setParameter("applyVars",&m_applyvars_str);
}

//////////////////////////////////////////////////////////////////////////////

BCFunction::~BCFunction()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCFunction::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< std::vector< std::string > >("Def","Definition of the functions (default <>)");
  options.addConfigOption< std::vector< std::string > >("applyVars","Names of variables to apply to (default <>)");
}

//////////////////////////////////////////////////////////////////////////////

void BCFunction::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  BC::configure(args);

  // check input correctness
  if (m_applyvars_str.size()!=m_function_str.size())
    err("\"applyVars\" and \"Def\" sizes should be the same");
}

//////////////////////////////////////////////////////////////////////////////

void BCFunction::parse()
{
  CFAUTOTRACE;

  // applyEqs
  m_applyvars.assign(m_applyvars_str.size(),-1);
  for (CFuint i=0; i<m_applyvars_str.size(); ++i) {
    m_applyvars[i] = getMethodData().getVariableIndex(m_applyvars_str[i]);
    if (m_applyvars[i]<0)
      err("variable \"" + m_applyvars_str[i] + "\" not found!");
  }

  // vectorial function
  m_function.setFunctions(m_function_str);
  m_function.setVariables(getMethodData().getNodalVariables());
  try {
    m_function.parse();
  }
  catch (ParserException& e) {
    ver("VectorialFunction parsing: " + std::string(e.what()));
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFunction::apply(const SafePtr< System > s)
{
  if (!m_function.isParsed())
    parse();

  const DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();
  const DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  RealVector vars(getMethodData().getNodalVariables().size());
  RealVector res(m_applyvars.size());
  const int iv   = s->iv;
  const int Nsys = s->Nsys;

  // for all the boundary regions, cycle the nodes
  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<trs.size(); ++i) {
    const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();
    for (std::vector< CFuint >::const_iterator
         n=nodes.begin(); n!=nodes.end(); ++n) {
      if (h_mn_priority[*n]<=m_bnpriority) {

        // evaluate VectorialFunction and set condition
        getMethodData().getNodalValues(*n,vars);
        m_function.evaluate(vars,res);
        for (CFuint i=0; i<m_applyvars.size(); ++i) {
          if (m_applyvars[i]>=iv && m_applyvars[i]<iv + Nsys)
            setDirichletCondition(*s,*n,m_applyvars[i]-iv,
              (*h_states[*n])[m_applyvars[i]] - res[i] );
        }

      }
    }
  }

  // apply to node index
  if (m_point_index>=0) {
    if (h_mn_priority[m_point_index]<=m_bnpriority) {

      // evaluate VectorialFunction and set condition
      getMethodData().getNodalValues(m_point_index,vars);
      m_function.evaluate(vars,res);
      for (CFuint i=0; i<m_applyvars.size(); ++i) {
        if (m_applyvars[i]>=iv && m_applyvars[i]<iv + Nsys)
          setDirichletCondition(*s,m_point_index,m_applyvars[i]-iv,
            (*h_states[m_point_index])[m_applyvars[i]] - res[i] );
      }

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFunction::setup()
{
  CFAUTOTRACE;
  BC::setup();

  if (!m_function.isParsed())
    parse();

  const DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();
  const DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  RealVector vars(getMethodData().getNodalVariables().size());
  RealVector res(m_applyvars.size());

  // for all the boundary regions, cycle the nodes
  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<trs.size(); ++i) {
    const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();
    for (std::vector< CFuint >::const_iterator
         n=nodes.begin(); n!=nodes.end(); ++n) {
      if (h_mn_priority[*n]<=m_bnpriority) {

        // evaluate VectorialFunction and set values
        getMethodData().getNodalValues(*n,vars);
        m_function.evaluate(vars,res);
        for (CFuint i=0; i<m_applyvars.size(); ++i)
          (*h_states[*n])[m_applyvars[i]] = res[i];

      }
    }
  }

  // apply to node index
  if (m_point_index>=0) {
    if (h_mn_priority[m_point_index]<=m_bnpriority) {

      // evaluate VectorialFunction and set values
      getMethodData().getNodalValues(m_point_index,vars);
      m_function.evaluate(vars,res);
      for (CFuint i=0; i<m_applyvars.size(); ++i)
        (*h_states[m_point_index])[m_applyvars[i]] = res[i];

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

