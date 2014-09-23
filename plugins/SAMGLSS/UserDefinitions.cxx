// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "samg.h"


namespace COOLFluiD {
  namespace SAMGLSS {

extern "C" {
#ifdef SAMG_PACKAGE_SAMGp

int __powi4i4(int i1, int i2)
{
  return i1^i2;
}

bool SAMG_USER_ALLOW_PARALLEL()
{
  return true;
}

SAMG_C_CALLCONV SAMG_USER_COO(int* i,int* ndim,double* x,double* y,double* z)
{
  //FIXME 2D and 3D coordinates
  *ndim=0;
}

#endif
#ifdef SAMG_PACKAGE_SAMG

bool SAMG_USER_ALLOW_PARALLEL()
{
  return false;
}

SAMG_C_CALLCONV SAMG_SET_REMOVE_REDUNDANCIES(int*)
{
}

SAMG_C_CALLCONV SAMG_SET_DO_DUMPRESTART(int*)
{
}

SAMG_C_CALLCONV SAMG_COMPTIME(float *t1, float *t2, float *t3)
{
  *t1 = *t3 - *t2;
}

SAMG_C_CALLCONV SAMGP(
  int *nnu, int *nna, int *nsys,
  int *ia, int *ja, double *a, double *f, double *u,
  int *iu, int *ndiu, int *ip, int *ndip,
  int *matrix, int *iscale,
  double *res_in, double *res_out, int *ncyc_done, int *ierr,
  int *nsolve, int *ifirst, double *eps, int *ncyc, int *iswtch,
  double *a_cmplx, double *g_cmplx, double *p_cmplx, double *w_avrge,
  double *chktol, int *idump, int *iout,
  int *nshalo, int *npsnd, int *iranksnd, int *ipts, int *isndlist,
  int *nrhalo, int *nprec, int *irankrec, int *iptr, int *ireclist,
  int *icomm)
{
  *ierr = 2;
}

#endif
} // extern "C"

  }  // namespace COOLFluiD
}  // namespace SAMGLSS

