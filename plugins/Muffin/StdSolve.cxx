
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Loop.hh"
#include "Muffin/StdSolve.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSolve,MuffinData,MuffinModule > cStdSolveProvider("StdSolve");

//////////////////////////////////////////////////////////////////////////////

void StdSolve::execute()
{
  // execute master loop
  std::vector< Common::SafePtr< Loop > >& vl = getMethodData().m_vcomm_loops;
  for (std::vector< Common::SafePtr< Loop > >::iterator l=vl.begin(); l!=vl.end(); ++l)
    if ((*l)->isMaster())
      (*l)->execute();

  // copy solution to rhs
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  for (CFuint n=0; n<h_states.size(); ++n)
    for (CFuint e=0; e<h_states[0]->size(); ++e)
      h_rhs(n,e,Neqns) = (*h_states[n])[e];
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

