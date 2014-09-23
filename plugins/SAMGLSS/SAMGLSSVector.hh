// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_SAMGLSS_SAMGLSSVector_hh
#define COOLFluiD_SAMGLSS_SAMGLSSVector_hh

#ifdef CF_HAVE_MPI
#include <mpi.h>
#endif
#include "Common/StringOps.hh"
#include "Framework/LSSVector.hh"
#include "MathTools/RealVector.hh"

namespace COOLFluiD {
  namespace SAMGLSS {

/// This class represents a proxy for the SAMGLSS vector
class SAMGLSSVector : public Framework::LSSVector {

public:

  /// Default constructor without arguments
  SAMGLSSVector() : Framework::LSSVector() {
    m_size = 0;
    m_v = CFNULL;
  }

  /// Destructor
  ~SAMGLSSVector() {
    destroy();
  }

  /// Create a vector
  void create(MPI_Comm comm, const CFint m, const CFint M, const char* name);

  /// Deallocate internal memory
  void destroy();

  /// Initialize a vector
  void initialize(MPI_Comm comm, const CFreal value) {
    setValue(value);
  }

  /// Start to assemble the vector
  void beginAssembly() {};

  /// Finish to assemble the vector
  void endAssembly() {};

  /// Print this vector
  void printToScreen() const;

  /// Print this vector to a file
  void printToFile(const char* fileName) const;


  /**
   * Set a value at the specified position in the vector
   * @param idx index in the vector
   * @param value value to be set
   */
  void setValue(const CFint idx, const CFreal value) {
    m_v[idx] = value;
  }

  /**
   * Set all the entries equal to the given value
   * Set all the entries equal to the given value
   * @param value value to be set everywhere
   */
  void setValue(const CFreal value) {
    for (CFuint i=0; i<m_size; ++i)
      m_v[i] = value;
  }

  /**
   * Set a list of values
   * Set a list of values
   * @param nbValues number of values to set
   * @param idx pointer to the index corresponding to the first value
   * @param values pointer to the first value
   */
  void setValues(
    const CFuint nbValues, const CFint* idx, const CFreal* values ) {
    for (CFuint i=0; i<nbValues; ++i)
      m_v[idx[i]] = values[i];
  }

  /**
   * Add a value in the vector at the given location
   * @param idx index in the vector
   * @param value value to be added
   */
  void addValue(const CFint idx, const CFreal value) {
    m_v[idx] += value;
  }

  /**
   * Add a list of values at the given locations
   * @param nbValues number of values to add
   * @param idx pointer to the index corresponding to the first value
   * @param values pointer to the first value
   */
  void addValues(
    const CFuint nbValues, const CFint* idx, const CFreal* values ) {
    for (CFuint i=0; i<nbValues; ++i)
      m_v[idx[i]] += values[i];
  }

  /// Get one value
  void getValue(const CFint idx, CFreal value) {
    value = m_v[idx];
  }

  /// Get a list of values
  void getValues(const CFuint m, const CFint* im, CFreal* values) {
    for (CFuint i=0; i<m; ++i)
      values[i] = m_v[im[i]];
  }


  /**
   * Copy the raw data of this vector to a given array
   * @param other  given array to which copy the content of m_v
   */
  void copy(CFreal *const other, const CFuint size) const {
    for (CFuint i=0; i<size; ++i)
      other[i] = m_v[i];
  }

  /**
   * Copy the raw data of this vector to a given array
   * @param other  given array to which copy the content of m_v
   */
  void copy(
    CFreal *const other, CFint *const localIDs, const CFuint size ) const {
    for (CFuint i=0; i<size; ++i)
      other[localIDs[i]] = m_v[i];
  }

  /// Gets the global size (corresponds to SAMG's NNU)
  CFuint getGlobalSize() const {
    return m_size_global;
  }


  // non-abstract member functions


  /// Get a single value
  double getValue(const CFint idx) {
    return m_v[idx];
  }

  /// Get internal array
  double* getArray() const {
    return &m_v[0];
  }

  /// Gets the local size of the vector
  CFuint getLocalSize() const {
    return m_size;
  }


private:

  /// Vector entries
  double* m_v;

  /// Size of local system (SAMG's nnu+nrhalo, where nrhalo=0 in serial)
  CFuint m_size;

  /// Size of global system (SAMG's NNU)
  CFuint m_size_global;

  /// Vector name
  std::string m_name;

  /// Copy constructor
  SAMGLSSVector(const SAMGLSSVector& other);

  /// Overloading of the assignment operator
  const SAMGLSSVector& operator= (const SAMGLSSVector& other);

}; // end of class SAMGLSSVector


  }  // namespace SAMGLSS
}  // namespace COOLFluiD

#endif // COOLFluiD_SAMGLSS_SAMGLSSVector_hh

