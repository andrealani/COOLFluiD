// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "Common/BadValueException.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/MethodCommandProvider.hh"

// compile without MPI, to check this is not a parallel simulation
#include "Common/PE.hh"

#include "Pardiso/PardisoConfig.hh"
#include "Pardiso/PardisoData.hh"
#include "Pardiso/PardisoModule.hh"

#ifdef CF_HAVE_PARDISO_UBASEL
#define PARDISO pardiso_
// U. Basel PARDISO function prototype is here because no header comes with
// distribution
extern "C" {
  void pardiso_( void*, int*, int*, int*, int*, int*, double*, int*, int*,
    int*, int*, int*, int*, double*, double*, int* );
}
#else
// MKL PARDISO function prototype and use mkl_set_dynamic/mkl_set_num_threads
// to control dynamic behaviour
#include "mkl_pardiso.h"
#include "mkl_service.h"
#endif

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Pardiso {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NullMethodCommand< PardisoData >, PardisoData,
  PardisoModule > nullPardisoComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void PardisoData::defineConfigOptions(Config::OptionList& options)
{
  // configuration parameters
  options.addConfigOption< bool >("OncePhase11","If phase 11 (analysis) is to happen only once (default false)");
  options.addConfigOption< bool >("OncePhase22","If phase 22 (numerical factorization) is to happen only once (default false)");
  options.addConfigOption< int >("MatrixType","Matrix type (1, 2, -2 or 11, default)");
  options.addConfigOption< int >("CGS","CGS LU-preconditioning control (0: direct LU preconditioner (default), other values (31,61,...) set CGS with LU preconditioner");
  options.addConfigOption< int >("OOC","In/out-of-core version selection (0: in-core (default), 1: automatic or 2: out-of-core)");
}

//////////////////////////////////////////////////////////////////////////////

PardisoData::PardisoData(Common::SafePtr< std::valarray< bool > > maskArray,
			 CFuint& nbSysEquations,
			 Common::SafePtr< Framework::Method > owner ) :
  LSSData(maskArray,nbSysEquations,owner),
  m_done_phase_11(false),
  m_done_phase_22(false)
{
  addConfigOptionsTo(this);

  if (Common::PE::GetPE().IsParallel()) {
    std::string e("PARDISO doesn't support MPI parallel solving");
    CFLog(ERROR,e << "\n");
    throw Common::NotImplementedException (FromHere(),e);
  }

  // set phases control and matrix type
  m_once_phase_11 = false;
  m_once_phase_22 = false;
  m_mtype = 11;
  setParameter("OncePhase11",&m_once_phase_11);
  setParameter("OncePhase22",&m_once_phase_22);
  setParameter("MatrixType",&m_mtype);

  // set control parameters (inc. CGS and OOC)
  for (int i=0; i<64; ++i)
    m_iparm[i] = 0;
  m_iparm[ 0] = 1;
  setParameter("CGS",&m_iparm[ 3]);
  setParameter("OOC",&m_iparm[59]);
}

//////////////////////////////////////////////////////////////////////////////

PardisoData::~PardisoData()
{
}

//////////////////////////////////////////////////////////////////////////////

void PardisoData::printToFile(const std::string& prefix, const std::string& suffix)
{
  CFAUTOTRACE;
  std::string matStr = prefix + "mat" + suffix;
  std::string rhsStr = prefix + "rhs" + suffix;
  std::string solStr = prefix + "sol" + suffix;

  m_mat.printToFile(matStr.c_str());
  m_rhs.printToFile(rhsStr.c_str());
  m_sol.printToFile(solStr.c_str());
}

//////////////////////////////////////////////////////////////////////////////

void PardisoData::configure(Config::ConfigArgs& args)
{
  LSSData::configure(args);

  // set phases control and matrix type
  //  1: real and structurally symmetric, supernode pivoting
  //  2: real and symmetric positive definite
  // -2: real and symmetric indefinite, diagonal or Bunch-Kaufman pivoting
  // 11: real and nonsymmetric, complete supernode pivoting
  m_once_phase_11 = (m_once_phase_22? true : m_once_phase_11);
  if (m_mtype!=1 && m_mtype!=2 && m_mtype!=-2 && m_mtype!=11) {
    std::string msg("PARDISO matrix type not recognized");
    CFLog(ERROR,msg << "\n");
    throw Common::BadValueException(FromHere(),msg);
  }

  // initialize internal memory pointer on first call to the solver only
  for (unsigned i=0; i<64; ++i)
    m_pt[i] = 0;

#ifndef CF_HAVE_PARDISO_UBASEL
  /*
   * WARNING. The current OOC version does not use threading, thus both the
   * parameter iparm(3) and the variable MKL_NUM_THREADS must be set to 1
   * when iparm(60) is equal to 1 or 2.
  */
  mkl_set_dynamic(0);                      // prevents setting 1 thread as default
  mkl_set_num_threads(m_iparm[59]? 1:0 );  // single/multi-threaded
  m_iparm[2] = (m_iparm[59]? 1:0);         // currently is not used
#endif
}

//////////////////////////////////////////////////////////////////////////////

void PardisoData::PardisoSolve(double *b, double *x)
{
  CFAUTOTRACE;

  // parameters
  int maxfct = 1;  // maximum number of numerical factorizations
  int mnum   = 1;  // which factorization to use
  int nrhs   = 1;  // number of right hand sides
  int error  = 0;  // initialize error flag
  int aperm  = 0;  // fill-in reducing ordering permutation vector (not used)
  int msglvl = (isOutput()? 1:0);  // print statistical information

  // system matrix
  int n     = m_mat.m_nnu;
  int *ia   = m_mat.getArrayIA();
  int *ja   = m_mat.getArrayJA();
  double *a = m_mat.getArrayA();

  // phase control
  std::vector< int > phases;
  if (!(m_once_phase_11 && m_done_phase_11))
    phases.push_back(11);
  if (!(m_once_phase_22 && m_done_phase_22))
    phases.push_back(22);
  phases.push_back(33);

  for (CFuint i=0; i<phases.size(); ++i) {
    int phase = phases[i];

    if (msglvl)
      CFLog(INFO,"PARDISO phase " << phase << "...\n");
    PARDISO( m_pt, &maxfct, &mnum, &m_mtype, &phase, &n, a, ia, ja, &aperm,
      &nrhs, m_iparm, &msglvl, b, x, &error );
    if (msglvl)
      CFLog(INFO,"PARDISO phase " << phase << ".\n");

    // control phases
    m_done_phase_11 = (phase==11? true:m_done_phase_11);
    m_done_phase_22 = (phase==22? true:m_done_phase_22);

    // output diagnostics
    if (phase==11 && msglvl) {
      CFLog(INFO,"PARDISO peak memory symbolic factorization:      " << m_iparm[14] << " kB\n");
      CFLog(INFO,"PARDISO permanent memory symbolic factorization: " << m_iparm[15] << " kB\n");
    }
    if (phase==22 && msglvl) {
      CFLog(INFO,"PARDISO memory numerical factorization and solution: " << m_iparm[16] << " kB\n");
      CFLog(INFO,"PARDISO peak memory solver consumption (total):      " << +m_iparm[14]+m_iparm[15]+m_iparm[16] << " kB\n");
    }
    if (phase==22 && m_iparm[3]) {
      if (m_iparm[19]>=0) {
        CFLog(INFO,"PARDISO CGS succeded with " << m_iparm[19] << " iterations\n");
      }
      else {
        CFLog(ERROR,"PARDISO CGS failed with error " << m_iparm[19] << "\n");
        // m_iparm[19] = - (cgs_it*10 + cgs_error)
        const int cgs_error = abs(m_iparm[19])%10;
        const int cgs_it    = abs(m_iparm[19])/10;
        std::string msg("PARDISO CGS error: "+Common::StringOps::to_str(cgs_error));
        switch (cgs_error) {
          case  1: { msg += ": too large fluctuations of the residuum"; } break;
          case  2: { msg += ": ||dxmax it cgs/2|| too large (slow convergence)"; } break;
          case  3: { msg += ": stopping criterion not reached at max it cgs"; } break;
          case  4: { msg += ": perturbed pivots caused iterative refinement"; } break;
          case  5: { msg += ": factorization is to fast for this matrix; it is better to use the factorization method with IPARM(4) = 0."; } break;
        }
        CFLog(ERROR,msg << "\n");
        CFLog(ERROR,"PARDISO CGS iterations: " << cgs_it << "\n");
      }
    }

    // output errors
    if (error) {
      std::string msg("PARDISO phase "+Common::StringOps::to_str(phase)+" error "+Common::StringOps::to_str(error));
      switch (error) {
        case  -1: { msg += ": input inconsistent"; } break;
        case  -2: { msg += ": not enough memory"; } break;
        case  -3: { msg += ": reordering problem"; } break;
        case  -4: { msg += ": zero pivot, numerical factorization or iterative refinement problem"; } break;
        case  -5: { msg += ": unclassified (internal) error"; } break;
        case  -6: { msg += ": preordering failed (matrix types 11, 13 only)"; } break;
        case  -7: { msg += ": diagonal matrix problem"; } break;
        case  -8: { msg += ": 32-bit integer overflow problem"; } break;
        case  -9: { msg += ": not enough memory for OOC"; } break;
        case -10: { msg += ": problems with opening OOC temporary files"; } break;
        case -11: { msg += ": read/write problems with the OOC data file"; } break;
      }
      CFLog(ERROR,msg << "\n");
      throw Common::BadValueException(FromHere(),msg);
      break;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PardisoData::PardisoDeallocate()
{
  CFAUTOTRACE;

  // parameters
  int maxfct = 1;  // maximum number of numerical factorizations
  int mnum   = 1;  // which factorization to use
  int nrhs   = 1;  // number of right hand sides
  int error  = 0;  // initialize error flag
  int aperm  = 0;  // fill-in reducing ordering permutation vector (not used)
  int msglvl = (isOutput()? 1:0);  // print statistical information
  double ddum = 0.;                // high precision dummy

  // system matrix
  int n     = (int) m_mat.m_nnu;
  int *ia   = m_mat.getArrayIA();
  int *ja   = m_mat.getArrayJA();
  int phase = -1;

  if (msglvl)
    CFLog(INFO,"PARDISO phase " << phase << "...\n");
  PARDISO( m_pt, &maxfct, &mnum, &m_mtype, &phase, &n, &ddum, ia, ja, &aperm,
    &nrhs, m_iparm, &msglvl, &ddum, &ddum, &error );
  if (msglvl)
    CFLog(INFO,"PARDISO phase " << phase << ".\n");

  if (error)
    CFLog(ERROR,"PARDISO phase " << phase << " error " << error << "\n");

  /// reset phases history
  m_done_phase_11 = false;
  m_done_phase_22 = false;
}

//////////////////////////////////////////////////////////////////////////////

unsigned PardisoData::digit(unsigned number, unsigned place, unsigned radix) const
{
  for (unsigned i=0; i<place; ++i)
    number /= radix;      // ignore preceding (small)
  return number % radix;  // and anteceding digits (big)
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Pardiso
} // namespace COOLFluiD

