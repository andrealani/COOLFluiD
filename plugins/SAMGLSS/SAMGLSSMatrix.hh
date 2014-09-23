// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_SAMGLSS_SAMGLSSMatrix_hh
#define COOLFluiD_SAMGLSS_SAMGLSSMatrix_hh

#ifdef CF_HAVE_MPI
#include <mpi.h>
#endif

#include "Common/StringOps.hh"
#include "Framework/LSSMatrix.hh"

namespace COOLFluiD {

  namespace Framework {  class BlockAccumulator;  }

  namespace SAMGLSS {

/// This class represents a proxy for the SAMGLSS matrix
class SAMGLSSMatrix : public Framework::LSSMatrix
{

public:

  /// Default constructor without arguments
  SAMGLSSMatrix() : Framework::LSSMatrix() {
    m_ia = CFNULL;
    m_ja = CFNULL;
    m_a  = CFNULL;
    m_nnu = 0;
    m_nna = 0;
  }

  /// Destructor
  ~SAMGLSSMatrix() {
    if (m_nnu) {
      delete[] m_ia;
      delete[] m_ja;
      delete[] m_a;
    }
    m_nnu = 0;
    m_nna = 0;
  }

  /// Create a sequential sparse matrix
  void createSeqAIJ(
    const CFint m, const CFint n,
    const CFint nz, const CFint* nnz,
    const char* name = CFNULL);

  /// Create a sequential block sparse matrix
  void createSeqBAIJ(
    const CFuint blockSize,
    const CFint m, const CFint n, const CFint nz, const CFint* nnz,
    const char* name = CFNULL);

  #ifdef CF_HAVE_MPI
  /// Create a parallel sparse matrix
  void createParAIJ(
    MPI_Comm comm,
    const CFint m, const CFint n, const CFint M, const CFint N,
    const CFint dnz, const int* dnnz, const CFint onz, const int* onnz,
    const char* name = CFNULL);
  #endif

  #ifdef CF_HAVE_MPI
  /// Create a parallel block sparse matrix
  void createParBAIJ(
    MPI_Comm comm, const CFuint blockSize,
    const CFint m, const CFint n, const CFint M, const CFint N,
    const CFint dnz, const int* dnnz, const CFint onz, const int* onnz,
    const char* name = CFNULL);
  #endif


  /// Start to assemble the matrix
  void beginAssembly(LSSMatrixAssemblyType assemblyType) {};

  /// Finish to assemble the matrix
  void endAssembly(LSSMatrixAssemblyType assemblyType) {
    if (!m_preconstruct)
      compressStructure();
  }

  /// Print this matrix
  void printToScreen() const;

  /// Print this matrix to a file
  void printToFile(const char* fileName) const;

  /// Set all the elements of the vector equal to the given value
  void setValue(const CFint im, const CFint in, const CFreal value);

  /// Set a list of values
  void setValues(
    const CFuint m, const CFint* im, const CFuint n, const CFint* in,
    const CFreal* values );

  /// Add an array of values in the vector
  void addValue(const CFint im, const CFint in, const CFreal value);

  /// Add a list of values
  void addValues(
    const CFuint m, const CFint* im, const CFuint n, const CFint* in,
    const CFreal* values );

  /// Get one value
  void getValue(
    const CFint im, const CFint in, CFreal& value );

  /// Get a list of values
  void getValues(
    const CFuint m, const CFint* im, const CFuint n, const CFint* in,
    CFreal* values );

  /// Set the diagonal
  void setDiagonal(Framework::LSSVector& diag);

  /// Add to the diagonal
  void addToDiagonal(Framework::LSSVector& diag);

  /// Reset to 0 all the non-zero elements of the matrix
  void resetToZeroEntries();

  /// Set a list of values
  void setValues(const Framework::BlockAccumulator& acc);

  /// Add a list of values
  void addValues(const Framework::BlockAccumulator& acc);

  /// Freeze the matrix structure concerning the non zero locations
  void freezeNonZeroStructure();


  // non-abstract member functions


  /// Create the matrix column structure (sets m_ja)
  void createJAStructure(std::vector< std::vector< CFuint > >& nz);

  /**
   * Get the matrix properties, which are: number of rows, entries and
   * unknowns and CSR IA, JA and A
   */
  void getMatrix( int& nnu, int& nna, int& nsys, int& matrix, int*& ia,
    int*& ja, double*& a ) const;

  /**
   * Manually set the matrix properties, which are: number of rows, entries
   * and unknowns, matrix type and CSR IA, JA and A
   */
  void setMatrix( int nnu, int nna, int nsys, int matrix, int* ia, int* ja,
    double* a );

  /// Get number of matrix entries (non-zeros)
  int getNbNonZeros() const {
    return m_nna;
  }

  /// Get number of variables (number of system rows)
  int getNbRows() const {
    return m_nnu;
  }

  /// Transform (row,column) indices pairs to CSR matrix indices
  void getCSRIndices( const CFuint m,const CFint* im,const CFuint n,
    const CFint* in,int*& ik );

  /// Transform BlockAccumulator indices to CSR matrix indices
  void getBlockAccumulatorCSRIndices(
    const Framework::BlockAccumulator& acc,int*& ik );

  /// Get one value, returning by value
  double getValue( const CFint im, const CFint in ) const;

  /// Get a diagonal value, returning by value
  double getDiagonalValue( const CFint im ) const {
    return m_a[m_ia[im]-1];
  }

  /// Handles situation for acessing non-allocated elements
  void notFoundIndex(const int r, const int c);

  /// Compress memory taken by JA
  void compressStructure();

  /// operator to get the matrix values by row and column indices
  /// @todo return a reference and remove non-abstract getValue(.,.)
  double operator() (const int r, const int c) {
    return getValue(r,c);
  }

  /// Manually set the matrix type
  void setMatrixType(int matrix) {
    m_matrix = matrix;
  }


private:

  /// Copy constructor
  SAMGLSSMatrix(const SAMGLSSMatrix& other);

  /// Overloading of the assignment operator
  const SAMGLSSMatrix& operator= (const SAMGLSSMatrix& other);

  /// Matrix name
  std::string m_name;

  /// Matrix rows indices (SAMG CSR ia)
  int* m_ia;

  /// Matrix column indices (SAMG CSR ja)
  int* m_ja;

  /// Matrix entries (SAMG CSR a)
  double* m_a;

  /// Matrix symmetry and zero row sum property (11, 12, 21 or 22)
  int m_matrix;

  /// Number of variables (system rows), same as updatable states * nbEqs
  int m_nnu;

  /// Number of system columns, same as local (upd. plus ghost) states * nbEqs
  int m_nnu_nrhalo;

  /// Number of matrix entries (non-zeros, generally)
  int m_nna;

  /// Number of unknowns, same as nbEqs
  int m_nsys;

  /// If the matrix structure is built in advance
  bool m_preconstruct;

}; // end of class SAMGLSSMatrix


  }  // namespace SAMGLSS
}  // namespace COOLFluiD

#endif // COOLFluiD_SAMGLSS_SAMGLSSMatrix_hh

