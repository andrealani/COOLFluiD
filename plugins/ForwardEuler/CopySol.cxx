// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/ForwardEuler.hh"

#include "ForwardEuler/CopySol.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< CopySol,FwdEulerData,ForwardEulerLib >
  copySolProvider("StdCopySol");

//////////////////////////////////////////////////////////////////////////////

void CopySol::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs     = socket_rhs.getDataHandle();

  const CFuint nbEqs    = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  for (CFuint i=0; i<nbStates; ++i) {
    cf_assert(states[i]->isValid());
    if (!states[i]->isParUpdatable())
      continue;

    for (CFuint j=0; j<nbEqs; ++j)
      (*states[i])[j] = rhs(i,j,nbEqs);
  }

  // computation of the norm of the dU (a.k.a rhs)
  CFreal value = 0.;
  const CFuint j = getMethodData().getVarID();
  for (CFuint i=0; i<nbStates; ++i) {
     const CFreal tmp = rhs(i,j,nbEqs);
     value += tmp*tmp;
  }
  value = std::max(value,MathTools::MathConsts::CFrealEps());
  getMethodData().setNorm( log10(sqrt(value)) );

}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > > CopySol::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler

} // namespace COOLFluiD

