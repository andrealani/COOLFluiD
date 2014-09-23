// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "SAMGLSS/SAMGLSSData.hh"
#include "SAMGLSS/SAMGLSSModule.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/NotImplementedException.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace SAMGLSS {

MethodCommandProvider< NullMethodCommand< SAMGLSSData >,SAMGLSSData,
  SAMGLSSModule >  nullSAMGLSSComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSData::defineConfigOptions(Config::OptionList& options)
{
  // configuration parameters
  options.addConfigOption< bool >(   "MatrixZeroRowSum",    "If matrix row entries sum to zero (default = false)." );
  options.addConfigOption< bool >(   "MatrixIsSymmetric",   "If matrix is symmetric (default = false)." );
  options.addConfigOption< bool >(   "UseInputChecking",    "If system should be checked and possibly fixed for library requirements compliance (default = true)." );
  options.addConfigOption< bool >(   "Preconstruct",        "If the matrix structure is to be built before assembly (default = true)." );
  options.addConfigOption< bool >(   "Debug",               "If debugging functions are to be called (default = false)." );
  options.addConfigOption< int >(      "OutputSetupLevel",    "Print output during the setup phase (idump parameter, default = 0)." );
  options.addConfigOption< int >(      "OutputSolutionLevel", "Print output during the setup phase (iout parameter, default = 2)." );
  options.addConfigOption< double >(   "RelativeTolerance",   "Relative tolerance for control of iterative solver convergence default = 1.0e-8)." );

  //FIXME these primary configuration parameters
  options.addConfigOption< int > ( "nsolve",    "nsolve (default = 2)." );
  options.addConfigOption< int > ( "ifirst",    "ifirst (default = 0)." );
  options.addConfigOption< int > ( "ncyc",      "ncyc (default = 11050)." );
  options.addConfigOption< int > ( "n_default", "n_default (default = 12)." );
  options.addConfigOption< int > ( "iswtch",    "iswtch (default = 5100+n_default)." );

  // primary parameters
  options.addConfigOption< CFreal >( "a_cmplx", "Memory upper limit for AMG operator complexity, tipically 1.5-3.0 (default = 2.2)." );
  options.addConfigOption< CFreal >( "g_cmplx", "Memory upper limit for AMG grid complexity, usually 2.0 (default = 1.7)." );
  options.addConfigOption< CFreal >( "p_cmplx", "Memory upper limit for AMG point complexity (irrelevant for scalar case), usually 2.0 (default = 2.0)." );
  options.addConfigOption< CFreal >( "w_avrge", "Memory upper limit for AMG average row length of interpolation, usually 3.0 (default = 2.4)." );

  // secondary parameters
  options.addConfigOption< bool >( "UseSecondaryParameters", "If secondary parameters should be manually set (default = false).");
  options.addConfigOption< CFint >(  "LEVELX", "Secondary parameter LEVELX (integer, default = 25)");
  options.addConfigOption< CFint >(  "NPTMN",  "Secondary parameter NPTMN (integer, default = 100)");
  options.addConfigOption< CFint >(  "NCG",    "Secondary parameter NCG (integer, default = 4)");
  options.addConfigOption< CFint >(  "NWT",    "Secondary parameter NWT (integer, default = 2)");
  options.addConfigOption< CFint >(  "NTR",    "Secondary parameter NTR (integer, default = 1)");
  options.addConfigOption< CFint >(  "NRD",    "Secondary parameter NRD (integer, default = 131)");
  options.addConfigOption< CFint >(  "NRU",    "Secondary parameter NRU (integer, default = 131)");
  options.addConfigOption< CFint >(  "NRC",    "Secondary parameter NRC (integer, default = 0)");
  options.addConfigOption< CFint >(  "NP_OPT", "Secondary parameter NP_OPT (integer, default = 0)");
  options.addConfigOption< CFreal >( "ECG",    "Secondary parameter ECG (real, default = 21.25)");
  options.addConfigOption< CFreal >( "EWT",    "Secondary parameter EWT (real, default = 0.20)");
  options.addConfigOption< CFreal >( "ETR",    "Secondary parameter ETR (real, default = 12.20)");
}

//////////////////////////////////////////////////////////////////////////////

SAMGLSSData::SAMGLSSData(Common::SafePtr< std::valarray< bool > > maskArray,
			 CFuint& nbSysEquations,
			 Common::SafePtr<Framework::Method> owner) :
  LSSData(maskArray, nbSysEquations, owner)
{
  addConfigOptionsTo(this);
  
  // solver configuration parameters
  m_sprimary.eps = 1.0e-8;
  m_sprimary.idump  = 0;
  m_sprimary.iout   = 2;
  setParameter( "RelativeTolerance",   &m_sprimary.eps );
  setParameter( "OutputSetupLevel",    &m_sprimary.idump );
  setParameter( "OutputSolutionLevel", &m_sprimary.iout );
  m_sprimary.a_cmplx = 2.2;
  m_sprimary.g_cmplx = 1.7;
  m_sprimary.p_cmplx = 2.0;
  m_sprimary.w_avrge = 2.4;
  setParameter( "a_cmplx", &m_sprimary.a_cmplx );
  setParameter( "g_cmplx", &m_sprimary.g_cmplx );
  setParameter( "p_cmplx", &m_sprimary.p_cmplx );
  setParameter( "w_avrge", &m_sprimary.w_avrge );

  //FIXME these primary configuration parameters
  // n_default configuration parameter
  // secondary parameters defaults: CURRENTLY AVAILABLE: 10-13, 15-18, 20-23,
  // 25-28. NOTE: the higher the respective SECOND digit, the more aggressive
  // the coarsening (lower memory at the expense of slower convergence).
  // secondary parameters have to be set if n_default=0.
  // iswtch configuration parameter:
  // memory de-allocation upon return; memory extension feature activated;
  // residuals measured in the L2-norm; secondary parameter default (setting
  // value of n_default)
  m_sprimary.nsolve = 2;      // solution strategy (6.4.1)
  m_sprimary.ifirst = 0;      // first approximation (0:u=u, 1:u=0, 2:u=1, 3:u=?)
  m_sprimary.ncyc   = 11050;  // V-cycle as CG pre-conditioner (at most 50 iterations)
  m_sprimary.n_default = 12;
  m_sprimary.iswtch = 5100+m_sprimary.n_default;
  setParameter( "nsolve",    &m_sprimary.nsolve );
  setParameter( "ifirst",    &m_sprimary.ifirst );
  setParameter( "ncyc",      &m_sprimary.ncyc );
  setParameter( "n_default", &m_sprimary.n_default);
  setParameter( "iswtch",    &m_sprimary.iswtch);

  // solver configuration parameters requiring processing
  bool matrix_zerorowsum = false;
  bool matrix_symmetric  = false;
  m_use_input_checking     = true;
  m_preconstruct           = true;
  m_debug                  = false;
  setParameter( "MatrixZeroRowSum",  &matrix_zerorowsum );
  setParameter( "MatrixIsSymmetric", &matrix_symmetric );
  setParameter( "UseInputChecking",  &m_use_input_checking );
  setParameter( "Preconstruct",      &m_preconstruct );
  setParameter( "Debug",             &m_debug );
  m_sprimary.matrix = (matrix_symmetric? 10:20)+(matrix_zerorowsum? 1:2);
  m_sprimary.chktol = (m_use_input_checking? 0.:-1.);  // 0.: logical check
  m_mat.setMatrixType(m_sprimary.matrix);

  // secondary parameters
  m_use_secondary_parameters = false;
  m_ssecondary.LEVELX =  25;
  m_ssecondary.NPTMN  = 100;
  m_ssecondary.NCG    =   4;
  m_ssecondary.NWT    =   2;
  m_ssecondary.NTR    =   1;
  m_ssecondary.NRD    = 131;
  m_ssecondary.NRU    = 131;
  m_ssecondary.NRC    =   0;
  m_ssecondary.NP_OPT =   0;
  m_ssecondary.ECG = 21.25;
  m_ssecondary.EWT =  0.20;
  m_ssecondary.ETR = 12.20;
  setParameter( "UseSecondaryParameters", &m_use_secondary_parameters );
  setParameter( "LEVELX", &m_ssecondary.LEVELX );
  setParameter( "NPTMN",  &m_ssecondary.NPTMN );
  setParameter( "NCG",    &m_ssecondary.NCG );
  setParameter( "NWT",    &m_ssecondary.NWT );
  setParameter( "NTR",    &m_ssecondary.NTR );
  setParameter( "NRD",    &m_ssecondary.NRD );
  setParameter( "NRU",    &m_ssecondary.NRU );
  setParameter( "NRC",    &m_ssecondary.NRC );
  setParameter( "NP_OPT", &m_ssecondary.NP_OPT );
  setParameter( "ECG",    &m_ssecondary.ECG );
  setParameter( "EWT",    &m_ssecondary.EWT );
  setParameter( "ETR",    &m_ssecondary.ETR );

  // set members
  m_counter = 0;

  m_sprimary.ia = CFNULL;
  m_sprimary.ja = CFNULL;
  m_sprimary.a  = CFNULL;
  m_sprimary.f  = CFNULL;
  m_sprimary.u  = CFNULL;
  m_sprimary.iu     = CFNULL;
  m_sprimary.ip     = CFNULL;
  m_sprimary.iscale = CFNULL;

  m_soutput.res_in  = 0.;
  m_soutput.res_out = 0.;
  m_soutput.ncyc_done = 0;
  m_soutput.ierr      = 0;
  m_soutput.told = 0.;
  m_soutput.tnew = 0.;
  m_soutput.tamg = 0.;

  m_sparallel.nrhalo = 0;
}

//////////////////////////////////////////////////////////////////////////////

SAMGLSSData::~SAMGLSSData()
{
  if (m_sprimary.iu!=CFNULL)      delete[] m_sprimary.iu;
  if (m_sprimary.ip!=CFNULL)      delete[] m_sprimary.ip;
  if (m_sprimary.iscale!=CFNULL)  delete[] m_sprimary.iscale;

  m_sprimary.iu     = CFNULL;
  m_sprimary.ip     = CFNULL;
  m_sprimary.iscale = CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSData::printToFile(const std::string& prefix, const std::string& suffix)
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

bool SAMGLSSData::isParallel() const
{
  bool b = Common::PE::GetPE().IsParallel();
  if (b && !SAMG_USER_ALLOW_PARALLEL()) {
    CFout << "SAMGLSSData: SAMG package doesn't support parallel solving\n";
    throw Common::NotImplementedException (FromHere(),
      "SAMGLSSData: SAMG package doesn't support parallel solving" );
  }
  return b;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSData::configure ( Config::ConfigArgs& args )
{
  LSSData::configure(args);
//  CFLog(NOTICE, "\n");
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSData::checkDataStructure()
{
  if (!m_use_input_checking)  return;

  CFLog(INFO, "SAMGLSS: checking system matrix and vectors...\n");
  m_sprimary_struct p = getParamsPrimary();

  // check and correct diagonal element positiveness. if negative is found,
  // reverse system row sign
  for (int k=0;k<p.nnu;++k)
    if (p.a[p.ia[k]-1]<0.) {
      for (int i=p.ja[k];i<p.ja[k+1]-1;++i)
        p.a[i] = -p.a[i];
      p.u[k] = -p.u[k];
      p.f[k] = -p.f[k];
    }

  // other tests (and fixes) can be added
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSData::getErrorDescription(int ierr,std::string& str)
{
  switch (ierr) {
  case   1:  str = "general error (unclassified)"; break;
  case   5:  str = "error in license checking"; break;
  case  10:  str = "insufficient dimensioning (happens only if memory extension is turned off)"; break;
  case  20:  str = "illegal input parameter"; break;
  case  30:  str = "undefined or missing input"; break;
  case  40:  str = "error in input arrays"; break;
  case  50:  str = "incorrect or inconsistent input"; break;
  case  51:  str = "re-start terminated: initial call with iswit=4 required"; break;
  case  52:  str = "re-start terminated: no AMG decomposition available"; break;
  case  53:  str = "re-start terminated: illegal change in parameters"; break;
  case  54:  str = "re-start terminated: illegal change in matrix"; break;
  case  60:  str = "memory management failed (including file I/O to scratch files)"; break;
  case  70:  str = "allocation and de-allocation errors"; break;
  case  71:  str = "unexpected allocation status of allocatable array"; break;
  case  72:  str = "unexpected assoziation status of pointer array"; break;
  case  80:  str = "requested AMG component no tinstalled in current release"; break;
  case  90:  str = "logfile exists but could not be opened"; break;
  case  91:  str = "logfile already connected to different unit"; break;
  case  92:  str = "logfile does not exist and could not be opened"; break;
  case  93:  str = "specified unit does not exist"; break;
  case  94:  str = "unit=5: specify a logfile or another unit number"; break;
  case 100:  str = "setup: general error (unclassified)"; break;
  case 101:  str = "setup: setup failed in CNTRL-routine (automatic setup mechanism)"; break;
  case 110:  str = "setup: error in defining strong connectivity"; break;
  case 120:  str = "setup: error in the splitting process"; break;
  case 130:  str = "setup: error in the coarsening process"; break;
  case 140:  str = "setup: error in defining interpolation"; break;
  case 150:  str = "setup: error in computing the coarse-level Galerkin operators"; break;
  case 160:  str = "setup: error in performing the setup optimization"; break;
  case 200:  str = "solution phase: general error (unclassified)"; break;
  case 210:  str = "solution phase: divergence of the method"; break;
  case 220:  str = "solution phase: error in relaxation (smoothing)"; break;
  case 230:  str = "solution phase: error in ILU (smoothing)"; break;
  case 240:  str = "solution phase: error in ILUT (smoothing)"; break;
  case 250:  str = "solution phase: error in inter-grid transfers"; break;
  case 260:  str = "solution phase: error in alternating \"Schwarz process\""; break;
  case 300:  str = "MPI-parallel: general error in communication"; break;
  case 310:  str = "MPI-parallel: illegal use of routine in parallel context"; break;
  case 320:  str = "MPI-parallel: illegal use of routine in sequential context"; break;
  case 398:  str = "MPI-parallel: a feature has been selected which is not yet supported by SAMGp"; break;
  case 399:  str = "MPI-parallel: a parameter has been accessed which is not yet supported by SAMGp"; break;
  case 800:  str = "auxiliary components: general error (unclassified)"; break;
  case 810:  str = "auxiliary components: error in ILU (one-level)"; break;
  case 820:  str = "auxiliary components: error in ILUT (one-level)"; break;
  case 830:  str = "auxiliary components: error in conjugate gradient (CG)"; break;
  case 831:  str = "auxiliary components: quasi residual check has been reached (CG)"; break;
  case 832:  str = "auxiliary components: residual stagnation check has been reached (CG)"; break;
  case 833:  str = "auxiliary components: limit of accelerator re-starts has been reached (CG)"; break;
  case 840:  str = "auxiliary components: error in BiCGstab"; break;
  case 841:  str = "auxiliary components: quasi residual check has been reached (BiCGstab)"; break;
  case 842:  str = "auxiliary components: residual stagnation check has been reached (BiCGstab)"; break;
  case 843:  str = "auxiliary components: limit of accelerator re-starts has been reached (BiCGstab)"; break;
  case 855:  str = "auxiliary components: error in GMRES"; break;
  case 851:  str = "auxiliary components: quasi residual check has been reached (GMRES)"; break;
  case 852:  str = "auxiliary components: residual stagnation check has been reached (GMRES)"; break;
  case 853:  str = "auxiliary components: limit of accelerator re-starts has been reached (GMRES)"; break;
  case 900:  str = "solution on coarsest level: general error (unclassified)"; break;
  case 910:  str = "solution on coarsest level: error in method #1 (iterative application of current smoother)"; break;
  case 920:  str = "solution on coarsest level: error in method #2 or #4 (preconditioned CG)"; break;
  case 930:  str = "solution on coarsest level: error in method #3 or #5 (preconditioned BiCGstab)"; break;
  case 960:  str = "solution on coarsest level: error in method #6 (full Gauss elimination)"; break;
  case 970:  str = "solution on coarsest level: error in method #7 (sparse Gauss elimination)"; break;
  case 980:  str = "solution on coarsest level: error in method #8 (least squares solver)"; break;
  case 990:  str = "solution on coarsest level: error in method #9 (user defined solver)"; break;

  // plugin-defined errors
  case   0:  str = "success"; break;
  case   2:  str = "function not available in library"; break;
  default:   str = "no description available"; break;
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

