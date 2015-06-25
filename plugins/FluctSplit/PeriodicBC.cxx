#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/PeriodicBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicBC,
                      FluctuationSplitData,
                      FluctSplitModule>
periodicBCProvider("PeriodicBC");

//////////////////////////////////////////////////////////////////////////////

void PeriodicBC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string > ("CoupledTrs","TRS to which this boundary is coupled.");
  options.addConfigOption< std::vector< std::string> >("Transform","Functions definind the transformation of coordinates");
  options.addConfigOption< CFreal >("Threshold","Threshold of distance that will consider two states matching after the transformation of coordinates");
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBC::PeriodicBC(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  m_match_states_idx(0)
{
  addConfigOptionsTo(this);

  m_transform_funcs = std::vector<std::string>();
  setParameter("Transform",&m_transform_funcs);

  m_coupled_trs = "";
  setParameter("CoupledTrs",&m_coupled_trs);

  m_threshold = 1e-12;
  setParameter("Threshold",&m_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBC::~PeriodicBC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PeriodicBC::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBC::setup()
{

  FluctuationSplitCom::setup();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  m_tcoord.resize(dim);
  m_delta.resize(dim);

  m_tmp_rhs.resize(nbEqs);

  // this will ensure trs exists
  Common::SafePtr<TopologicalRegionSet> coupled_trs =
    MeshDataStack::getActive()->getTrs(m_coupled_trs);

  // make sure that this boundary is only applied to one trs
  const std::vector<std::string>& trs_names = this->getTrsNames();
  if (trs_names.size() != 1)
  {
    std::string msg = "Periodic boundary condition" + getName() + " is not applying to only one TRS";
    throw BadValueException (FromHere(),msg);
  }

  // test if both TRS's have the same number of states
  Common::SafePtr<TopologicalRegionSet> applied_trs =
    MeshDataStack::getActive()->getTrs(trs_names[0]);

  Common::SafePtr< vector<CFuint> > const applied_trs_statesIdx = applied_trs->getStatesInTrs();
  Common::SafePtr< vector<CFuint> > const coupled_trs_statesIdx = coupled_trs->getStatesInTrs();

  const CFuint nb_states_applied_trs = applied_trs_statesIdx->size();
  const CFuint nb_states_coupled_trs = coupled_trs_statesIdx->size();

  if ( nb_states_applied_trs != nb_states_coupled_trs )
  {
    std::string msg = "Periodic boundary condition " + getName() +
                   " is applied to TRS "  + applied_trs->getName() +
                   " and coupled to TRS " + coupled_trs->getName() +
                   " which have different number of states " +
                   StringOps::to_str(nb_states_applied_trs) + " and " +
                   StringOps::to_str(nb_states_coupled_trs) + "\n";

    throw BadValueException (FromHere(),msg);
  }

  // build the linking vector of nodes
  m_match_states_idx.resize(nb_states_applied_trs);

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFreal sq_thres = m_threshold * m_threshold;

  for (CFuint is = 0; is < nb_states_applied_trs; ++is)
  {
    bool matched = false;

    const CFuint applied_sid = (*applied_trs_statesIdx)[is];
    State *const applied_state = states[applied_sid];
    Node& applied_coord = applied_state->getCoordinates();

//     CF_DEBUG_OBJ(applied_coord);

    m_vFunction.evaluate(applied_coord,m_tcoord);

    for (CFuint cs = 0; cs < nb_states_coupled_trs; ++cs)
    {
      const CFuint coupled_sid = (*coupled_trs_statesIdx)[cs];
      State *const coupled_state = states[coupled_sid];
      Node& coupled_coord = coupled_state->getCoordinates();

      m_delta = coupled_coord - m_tcoord;

      if ( m_delta.sqrNorm() < sq_thres ) // comparing the squares is faster
      {
        m_match_states_idx[is] = coupled_sid;
        matched = true;
        break; // go to next applied state
      }
    }

    if (!matched)
    {
      std::cout<<"PeriodicBC was not able to match all states between applied and coupled TRS"<<std::endl;
      throw BadValueException (FromHere(),"PeriodicBC was not able to match all states between applied and coupled TRS");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBC::executeOnTrs()
{
// unused //  DataHandle < Framework::State*, Framework::GLOBAL > states      = socket_states.getDataHandle();
// unused //  DataHandle<bool>    isUpdated   = socket_isUpdated.getDataHandle();
  DataHandle<CFreal>  rhs         = socket_rhs.getDataHandle();
  DataHandle<CFreal>  updateCoeff = socket_updateCoeff.getDataHandle();
 // get the data handle for the is Updated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  Common::SafePtr< vector<CFuint> > const applied_trs_statesIdx = getCurrentTRS()->getStatesInTrs();

  const CFuint nbEqs = m_tmp_rhs.size();

  for (CFuint is = 0; is < applied_trs_statesIdx->size(); ++is)
  {
    const CFuint applied_sid = (*applied_trs_statesIdx)[is];
   
// unused //    State *const applied_state = states[applied_sid];

    const CFuint coupled_sid = m_match_states_idx[is];
// unused //    State *const coupled_state = states[coupled_sid];

    // residuals
   

    if ((!isUpdated[applied_sid]) && (!isUpdated[coupled_sid])){

    for (CFuint j = 0; j < nbEqs; ++j)
        m_tmp_rhs[j] =  rhs(applied_sid,j,nbEqs);

    for (CFuint j = 0; j < nbEqs; ++j)
        m_tmp_rhs[j] += rhs(coupled_sid,j,nbEqs);

    for (CFuint j = 0; j < nbEqs; ++j)
    {
      rhs(coupled_sid,j,nbEqs) = m_tmp_rhs[j];
      rhs(applied_sid,j,nbEqs) = m_tmp_rhs[j];
    }

    // update coefficient

    const CFreal tmp_upcoeff = updateCoeff[applied_sid] + updateCoeff[coupled_sid];
    updateCoeff[applied_sid] = tmp_upcoeff;
    updateCoeff[coupled_sid] = tmp_upcoeff;
    isUpdated[applied_sid] = true;
    isUpdated[coupled_sid] = true;


  }

  }
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBC::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const CFuint dim = physModel->getDim();

  if ( m_transform_funcs.size() != dim )
  {
    throw BadValueException (FromHere(),"PeriodicBC : number of functions is not equal to dimensionality of the problem");
  }

  m_vars.resize(dim);
  if (dim >= DIM_1D) { m_vars[XX] = "x"; }
  if (dim >= DIM_2D) { m_vars[YY] = "y"; }
  if (dim >= DIM_3D) { m_vars[ZZ] = "z"; }

  // try to parse the transformation functions
  m_vFunction.setFunctions(m_transform_funcs);
  m_vFunction.setVariables(m_vars);
  try
  {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e)
  {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
