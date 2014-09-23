// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_SAMGLSS_SAMGLSSData_hh
#define COOLFluiD_SAMGLSS_SAMGLSSData_hh

#include "Framework/LSSData.hh"
#include "SAMGLSS/UserDefinitions.hh"
#include "SAMGLSS/SAMGLSSMatrix.hh"
#include "SAMGLSS/SAMGLSSVector.hh"

namespace COOLFluiD {
  namespace SAMGLSS {

//////////////////////////////////////////////////////////////////////////////

/// Structure holding the solver's primary parameters
struct m_sprimary_struct {
  int nnu,nna,nsys,*ia,*ja;
  double *a,*f,*u;
  int ndiu,ndip,*iu,*ip,*iscale;
  double eps,chktol,a_cmplx,g_cmplx,p_cmplx,w_avrge;
  int matrix,nsolve,ifirst,ncyc,n_default,iswtch,idump,iout;
};

/// Structure holding the solver's secondary parameters
struct m_ssecondary_struct {
  int LEVELX,NPTMN,NCG,NWT,NTR,NRD,NRU,NRC,NP_OPT;
  double ECG,EWT,ETR;
};

/// Structure holding the solver's parallel parameters
struct m_sparallel_struct {
  int icomm,
    nshalo,npsnd,*iranksnd,*ipts,*isndlist,
    nrhalo,nprec,*irankrec,*iptr,*ireclist;  // SAMGp
//  int *local2global,*partition,nnuglobal;    // extra SAMGpp
//  int ndlocal2global,ivar_min,ivar_max;      // extra SAMGp_pcrs
};

/// Structure holding the solver's output and timing parameters (6.9)
struct m_soutput_struct {
  double res_in,res_out;
  int ncyc_done,ierr;
  float told,tnew,tamg;
};

//////////////////////////////////////////////////////////////////////////////

/// This is a data object accessed by SAMGLSSCom's that compose the SAMGLSS
class SAMGLSSData : public Framework::LSSData
{

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the options
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  SAMGLSSData(Common::SafePtr< std::valarray< bool > > maskArray,
	      CFuint& nbSysEquations,
	      Common::SafePtr< Framework::Method > owner);
  
  /// Destructor
  ~SAMGLSSData();

  /// Gets if this is a parallel run, cf_asserting the library can solve
  bool isParallel() const;

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );


  /// Gets the SAMGLSS matrix
  SAMGLSSMatrix& getMatrix() {
    return m_mat;
  }

  /// Gets the SAMGLSS solution vector
  SAMGLSSVector& getSolVector() {
    return m_sol;
  }

  /// Gets the SAMGLSS RHS vector
  SAMGLSSVector& getRhsVector() {
    return m_rhs;
  }

  /// Set the SAMG interface to use
  void setInterface(const std::string& interface) {
    m_interface = interface;
  }

  /// Get the SAMG interface to use
  std::string getInterface() const {
    return m_interface;
  }

  /// Gets solver's primary parameters set
  m_sprimary_struct& getParamsPrimary() {
    m_mat.getMatrix(
      m_sprimary.nnu,m_sprimary.nna,m_sprimary.nsys,m_sprimary.matrix,
      m_sprimary.ia,m_sprimary.ja,m_sprimary.a );
    m_sprimary.f = m_rhs.getArray();
    m_sprimary.u = m_sol.getArray();
    return m_sprimary;
  }

  /// Gets solver's secondary parameters sst
  m_ssecondary_struct& getParamsSecondary() {
    return m_ssecondary;
  }

  /// Get solver's parallel parameters
  m_sparallel_struct& getParamsParallel() {
    return m_sparallel;
  }

  /// Get solver's output parameters
  m_soutput_struct& getParamsOutput() {
    return m_soutput;
  }

  /// Prints the Linear System to a file
  void printToFile(const std::string& prefix, const std::string& suffix);

  /// Gets the Class name
  static std::string getClassName() {
    return "SAMGLSS";
  }

  /// Set number of states present in the simulation (in a serial environment
  /// all the parameters have the same value)
  // (non-abstract member function)
  void setNbStates(CFuint nnb_states, CFuint nnb_local_states, CFuint
    nnb_global_states) {
    m_nb_states        = nnb_states;
    m_nb_local_states  = nnb_local_states;
    m_nb_global_states = nnb_global_states;
  }

  /// Get number of states present in the simulation (in a serial environment
  /// all the parameters have the same value)
  // (non-abstract member function)
  void getNbStates(CFuint& nnb_states, CFuint& nnb_local_states,
    CFuint& nnb_global_states) const {
    nnb_states        = m_nb_states;
    nnb_local_states  = m_nb_local_states;
    nnb_global_states = m_nb_global_states;
  }

  /// Returns if the matrix structure is to be built before assembly
  // (non-abstract member function)
  bool isPreconstruct() {
    return m_preconstruct;
  }

  /// Returns if debugging functions are to be called
  // (non-abstract member function)
  bool isDebug() {
    return m_debug;
  }

  /// Returns if secondary parameters are to be set manually (not-default
  /// parameters set)
  // (non-abstract member function)
  bool useSecondaryParameters() {
    return !m_sprimary.n_default;
  }

  /// Get and then increment counter of number of calls to solver
  // (non-abstract member function)
  CFuint getCounter() {
    return m_counter++;
  }

  /// Check data structure and attempt to fix errors (use just before solving)
  // (non-abstract member function)
  void checkDataStructure();

  /// Get error description based on SAMG return code
  // (non-abstract member function)
  void getErrorDescription(int ierr,std::string& str);


private:  // data

  /// System matrix
  SAMGLSSMatrix m_mat;

  /// System solution vector
  SAMGLSSVector m_sol;

  /// System right hand side (RHS) vector
  SAMGLSSVector m_rhs;

  /// Check system for library requirements compliance
  bool m_use_input_checking;

  /// SAMGp interface string
  std::string m_interface;

  /// If secondary parameters should be manually set
  bool m_use_secondary_parameters;

  /// Solver's primary parameters set
  m_sprimary_struct m_sprimary;

  /// Solver's secondary parameters set
  m_ssecondary_struct m_ssecondary;

  /// Solver's parallel parameters set
  m_sparallel_struct m_sparallel;

  /// Solver's parallel parameters set
  m_soutput_struct m_soutput;

  /// Number of states local to the computational node including non-updatable
  CFuint m_nb_states;

  /// Number of states local and updatable to the computational node (used
  /// only in parallel environment)
  CFuint m_nb_local_states;

  /// Number of states of the complete parallel system (used only in parallel
  /// environment)
  CFuint m_nb_global_states;

  /// If the matrix structure is to be built before assembly
  bool m_preconstruct;

  /// If debugging functions are to be called
  bool m_debug;

  /// Number of calls to solver
  CFuint m_counter;

}; // end of class SAMGLSSData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for SAMGLSS
typedef Framework::MethodCommand< SAMGLSSData > SAMGLSSCom;

/// Definition of a command provider for SAMGLSS
typedef Framework::MethodCommand< SAMGLSSData >::PROVIDER SAMGLSSComProvider;

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

#endif // COOLFluiD_SAMGLSS_SAMGLSSData_hh

