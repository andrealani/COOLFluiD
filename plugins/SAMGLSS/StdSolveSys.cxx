// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "samg.h"
#include "SAMGLSS/StdSolveSys.hh"
#include "SAMGLSS/SAMGLSSModule.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace SAMGLSS {

MethodCommandProvider< StdSolveSys,SAMGLSSData,SAMGLSSModule >
  stdSolveSysProvider("StdSolveSys");

//////////////////////////////////////////////////////////////////////////////

void StdSolveSys::execute()
{
  CFAUTOTRACE;

  // fix data structure according to library requirements
  getMethodData().checkDataStructure();

  // necessary before consecutive calls
  if (getMethodData().getCounter())
    SAMG_RESET_HIDDEN();

  m_sprimary_struct&   p = getMethodData().getParamsPrimary();
  m_sparallel_struct&  y = getMethodData().getParamsParallel();
  m_soutput_struct&    o = getMethodData().getParamsOutput();

  DataHandle<CFreal > h_rhs = socket_rhs.getDataHandle();
  LSSIdxMapping& L2S = getMethodData().getLocalToGlobalMapping();

  // copy h_rhs into samg rhs vector (p.f)
  const int nbRows = (p.nnu+y.nrhalo)/p.nsys;
  for (int iL=0;iL<nbRows;++iL) {            // local index
    const int iS = L2S.getColID(iL,p.nsys);  // SAMG index
    if (iS<p.nnu)
      for (int iEq=0;iEq<p.nsys;++iEq)
        p.f[iS+iEq] = h_rhs(iL,iEq,p.nsys);
  }

  std::string interface = getMethodData().getInterface();
  CFLog(INFO, "SAMGLSS: calling " << interface << "...\n");
  SAMG_CTIME(&o.told);
  if (interface=="SAMG") {
    SAMG(
      &p.nnu,&p.nna,&p.nsys,
      p.ia,p.ja,p.a,p.f,p.u,
      p.iu,&p.ndiu,p.ip,&p.ndip,&p.matrix,p.iscale,
      &o.res_in,&o.res_out,&o.ncyc_done,&o.ierr,
      &p.nsolve,&p.ifirst,&p.eps,&p.ncyc,&p.iswtch,
      &p.a_cmplx,&p.g_cmplx,&p.p_cmplx,&p.w_avrge,
      &p.chktol,&p.idump,&p.iout );
  } else if (interface=="SAMGp") {
    /* int iunformatted=0; SAMGP_C_DUMPRESTART( &iunformatted, */
    SAMGP(
      &p.nnu,&p.nna,&p.nsys,
      p.ia,p.ja,p.a,p.f,p.u,
      p.iu,&p.ndiu,p.ip,&p.ndip,&p.matrix,p.iscale,
      &o.res_in,&o.res_out,&o.ncyc_done,&o.ierr,
      &p.nsolve,&p.ifirst,&p.eps,&p.ncyc,&p.iswtch,
      &p.a_cmplx,&p.g_cmplx,&p.p_cmplx,&p.w_avrge,
      &p.chktol,&p.idump,&p.iout,
      &y.nshalo,&y.npsnd,y.iranksnd,y.ipts,y.isndlist,
      &y.nrhalo,&y.nprec,y.irankrec,y.iptr,y.ireclist,
      &y.icomm );
  }
  SAMG_CTIME(&o.tnew);
  SAMG_COMPTIME(&o.tamg,&o.told,&o.tnew);

  // copy samg solution vector (p.u) to h_rhs
  for (int iL=0;iL<nbRows;++iL) {
    const int iS = L2S.getColID(iL,p.nsys);
    if (iS<p.nnu)
      for (int iEq=0;iEq<p.nsys;++iEq)
        h_rhs(iL,iEq,p.nsys) = p.u[iS+iEq];
  }

  // output
  std::string err;
  getMethodData().getErrorDescription(std::abs(o.ierr),err);
  CFLog(INFO, "SAMGLSS: terminated with " << (o.ierr>0? "error ":"warning ")
    << o.ierr << ": " << err << "\n");
  CFLog(INFO, "SAMGLSS: converged at iteration: " << o.ncyc_done << "\n");
  CFLog(INFO, "SAMGLSS: time: " << o.tamg << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

