
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/SystemMITReM.hh"
#include "Muffin/BCBulk.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCBulk,MuffinData,MuffinModule > cBCBulkProvider("Bulk");

//////////////////////////////////////////////////////////////////////////////

BCBulk::BCBulk(const std::string& name) :
    BC(name),
    m_concentrations_recalculated(false)
{
  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCBulk");
  addConfigOptionsTo(this);

  m_factor = 1.;
  m_concentrations.clear();
  setParameter("Factor",&m_factor);
  setParameter("Concentrations",&m_concentrations);
}

//////////////////////////////////////////////////////////////////////////////

void BCBulk::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< double >("Factor","Factor to multiply bulk concentrations with (default 1.)");
  options.addConfigOption< std::vector< double > >("Concentrations","Concentrations to apply (overrides \"Factor\", default 1.)");
}

//////////////////////////////////////////////////////////////////////////////

void BCBulk::applyOnSystemMITReM(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  SystemMITReM& sys = dynamic_cast< SystemMITReM& >(*s);
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  // adjust concentrations to be electrically neutral (done here because
  // library is needed)
  if (!m_concentrations_recalculated) {

    if (!m_concentrations.size()) {
      m_concentrations = sys.m_bulk;
      for (int i=0; i<sys.Nions; ++i)
        m_concentrations[i] *= m_factor;
    }
    else {
      if ((int) m_concentrations.size()!=sys.Nions)
        err("\"Concentrations\" size must be number of ions");
      sys.m_mitrem->calcEquilibrium(m_concentrations);
    }

    log("set bulk concentrations...");
    for (CFuint i=0; i<sys.m_mitrem->getNIons(); ++i) {
      log(sys.m_mitrem->getIonLabel(i) + ": " + StringOps::to_str(m_concentrations[i]));
    }
    ver("set bulk concentrations.");
    m_concentrations_recalculated = true;
  }

  /*
   * new method: apply Dirichlet condition to all equations but the last,
   * which is the charge flux conservation. works beautifully.
   *
   * (older: the following doesn't exactly apply, it's here for reference)
   *
   * this applies the equivalent to a Dirichlet boundary condition for the
   * (bulk) ion species concentration, which should be in equilibrium as
   * calculated by the library. it's not possible to enforce identity to each
   * mass trasport equation because the remaining equation (electrostatics)
   * then is a linear combination of the other (for electroneutrality; for
   * Poisson this is 'numerically' true.) this makes the system singular (or
   * effectivelly badly conditioned as an iterative solver would see it).
   *
   * a way to get efficiently (but not exactly correct), is to set the mass
   * transport equations minus one, and the electrostatics equation, as a
   * Dirichlet condition for the species concentrations. this leaves the
   * remaining mass transport equation to couple the problem but it holds well
   * better then electrostatics (see formulation.)
   *
   * to implement, there are two remarks: since there is the correct (bulk)
   * value in the solution vector and this is a Newton method (we're solving
   * for the variation of the solution, this implies that the imposition is of
   * the form [0 0 0 0 1 0 0] [x] = [0].
   *
   * secondly, since we have to swap the system equations order to avoid a
   * zero in the diagonal -- that is, moving the electrostatics equation out
   * of the last system row, it suffices to apply the condition to the first
   * (number of mass transport equations) to make sure the electrostatics
   * equation is overwritten. also, it doesn't matter which variable the
   * equation will set (the column where the '1' sits) because the correct
   * solution values are already there anyway, so we might just as well put
   * them diagonally to help the solver.
   */

  // for all the boundary regions, cycle the nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {

    // for mass balance and electroneutrality equations, charge flux not set
    for (int i=0; i<sys.Nions; ++i)
      setDirichletCondition( sys,*n,i,
        (*h_states[*n])[sys.iv+i] - m_concentrations[i] );

  }

  // apply to node index
  if (m_point_index>=0) {
    for (int i=0; i<sys.Nions; ++i)
      setDirichletCondition( sys,m_point_index,i,
        (*h_states[m_point_index])[sys.iv+i] - m_concentrations[i] );
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

