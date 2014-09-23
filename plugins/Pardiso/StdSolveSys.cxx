// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

#include "Pardiso/StdSolveSys.hh"
#include "Pardiso/PardisoModule.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Pardiso {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSolveSys,PardisoData,PardisoModule >
  stdSolveSysProvider("StdSolveSys");

//////////////////////////////////////////////////////////////////////////////

void StdSolveSys::execute()
{
  PardisoData& d = getMethodData();

  const int Nstate = (int) d.getNbStates();     // nb. states
  const int Nsys   = (int) d.getNbEquations();  // nb. unknowns

  // create equation index mapping (solver to absolute)
  std::vector< int > equationid(Nsys);
  Common::SafePtr< std::valarray< bool > > maskarray = d.getMaskArray();
  const CFuint Neqns = maskarray->size();
  int e = 0;
  for (CFuint i=0; i<Neqns; ++i)
    if ((*maskarray)[i])
      equationid[e++] = i;
  cf_assert(e==Nsys);

  // system vectors
  double *b = d.getRhsVector().getArray();
  double *x = d.getSolVector().getArray();

  // copy h_rhs into right-hand side vector
  DataHandle<CFreal > h_rhs = socket_rhs.getDataHandle();
  for (int n=0; n<Nstate; ++n)
    for (int e=0; e<Nsys; ++e)
      b[n*Nsys+e] = h_rhs(n,equationid[e],Neqns);

  // call solver
  d.PardisoSolve(b,x);

  // copy solution to right-hand side vector
  for (int n=0; n<Nstate; ++n)
    for (int e=0; e<Nsys; ++e)
      h_rhs(n,equationid[e],Neqns) = x[n*Nsys+e];
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Pardiso
}  // namespace COOLFluiD

