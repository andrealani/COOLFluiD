// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Pardiso_PardisoData_hh
#define COOLFluiD_Pardiso_PardisoData_hh

#include "Framework/LSSData.hh"
#include "Pardiso/PardisoMatrix.hh"
#include "Pardiso/PardisoVector.hh"

namespace COOLFluiD {
  namespace Pardiso {

//////////////////////////////////////////////////////////////////////////////

/// This is a data object accessed by PardisoComs
class PardisoData : public Framework::LSSData
{

 public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the options
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  PardisoData(Common::SafePtr< std::valarray< bool > > maskArray,
    CFuint& nbSysEquations,
    Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~PardisoData();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Get the PardisoMatrix
  PardisoMatrix& getMatrix() {
    return m_mat;
  }

  /// Get the PardisoVector for the solution
  PardisoVector& getSolVector() {
    return m_sol;
  }

  /// Get the PardisoVector for the RHS
  PardisoVector& getRhsVector() {
    return m_rhs;
  }

  /// Prints the Linear System to a file
  void printToFile(const std::string& prefix, const std::string& suffix);

  /// Get the class name
  static std::string getClassName() {
    return "Pardiso";
  }

  /// Set number of states
  // (non-abstract member function)
  void setNbStates(CFuint nnb_states) {
    m_nb_states = nnb_states;
  }

  /// Get number of states
  // (non-abstract member function)
  CFuint getNbStates() const {
    return m_nb_states;
  }

  /// Set number of variables per state
  // (non-abstract member function)
  void setNbEquations(CFuint nnb_equations) {
    m_nb_equations = nnb_equations;
  }

  /// Get number of variables per state
  // (non-abstract member function)
  CFuint getNbEquations() const {
    return m_nb_equations;
  }

  /// Call PARDISO solver taking right-hand side and solution vectors
  // (non-abstract member function)
  void PardisoSolve(double *b, double *x);

  /// Call PARDISO for deallocation
  // (non-abstract member function)
  void PardisoDeallocate();


 private:  // methods

  /// Extract a digit from a positive integer, according to a certain radix
  unsigned digit(unsigned number, unsigned place, unsigned radix=10) const;


 private:  // data

  /// System matrix
  PardisoMatrix m_mat;

  /// System solution vector
  PardisoVector m_sol;

  /// System right hand side (RHS) vector
  PardisoVector m_rhs;

  /// If phase 11 (analysis) is to happen only once (configurable)
  bool m_once_phase_11;

  /// If phase 22 (numerical factorization) is to happen only once (configurable)
  bool m_once_phase_22;

  /// If phase 11 (analysis) was done before
  bool m_done_phase_11;

  /// If phase 22 (numerical factorization) was done before
  bool m_done_phase_22;

  /// Matrix type (configurable)
  int m_mtype;

  /// Number of states
  CFuint m_nb_states;

  /// Number of variables per state
  CFuint m_nb_equations;

  /// Internal memory pointer
  /// - 32-bit: int pt[64]
  /// - 64-bit: long int pt[64]
  /// (void *pt[64] should be OK on both architectures)
  void *m_pt[64];

  /// Control parameters
  int m_iparm[64];

}; // end of class PardisoData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for Pardiso
typedef Framework::MethodCommand< PardisoData > PardisoCom;

/// Definition of a command provider for Pardiso
typedef Framework::MethodCommand< PardisoData >::PROVIDER PardisoComProvider;

//////////////////////////////////////////////////////////////////////////////

  } // namespace Pardiso
} // namespace COOLFluiD

#endif // COOLFluiD_Pardiso_PardisoData_hh

