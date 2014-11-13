// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LSSMatrix_hh
#define COOLFluiD_Framework_LSSMatrix_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonCopyable.hh"

#include "Framework/Framework.hh"

#ifdef CF_HAVE_MPI
#  include <mpi.h>
#endif // CF_HAVE_MPI

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class BlockAccumulator;
    class LSSVector;

//////////////////////////////////////////////////////////////////////////////

/// This class represents an abstract Linear System Solver Matrix
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API LSSMatrix : public Common::NonCopyable<LSSMatrix> {
public:

  enum LSSMatrixAssemblyType {FLUSH_ASSEMBLY, FINAL_ASSEMBLY};

  /// Default constructor without arguments.
  LSSMatrix();

  /// Destructor.
  virtual ~LSSMatrix();
 
  /// use GPU support
  void setGPU(bool useGPU) {m_useGPU = useGPU;}
  
  /// Create a sequential sparse matrix
  virtual void createSeqAIJ(const CFint m,
                            const CFint n,
                            const CFint nz,
                            const CFint* nnz,
                            const char* name = CFNULL) = 0;

  /// Create a sequential block sparse matrix
  virtual void createSeqBAIJ(const CFuint blockSize,
                             const CFint m,
                             const CFint n,
                             const CFint nz,
                             const CFint* nnz,
                             const char* name = CFNULL) = 0;

#ifdef CF_HAVE_MPI

  /// Create a parallel sparse matrix
  virtual void createParAIJ(MPI_Comm comm,
                            const CFint m,
                            const CFint n,
                            const CFint M,
                            const CFint N,
                            const CFint dnz,
                            const CFint* dnnz,
                            const CFint onz,
                            const CFint* onnz,
                            const char* name = CFNULL) = 0;
  
  /// Create a parallel block sparse matrix
  virtual void createParBAIJ(MPI_Comm comm,
                             const CFuint blockSize,
                             const CFint m,
                             const CFint n,
                             const CFint M,
                             const CFint N,
                             const CFint dnz,
                             const CFint* dnnz,
                             const CFint onz,
                             const CFint* onnz,
                             const char* name = CFNULL) = 0;
#endif // CF_HAVE_MPI

  /// Perform an intermediate assembly when switching
  /// between ADD and SET values or viceversa
  void flushAssembly();

  /// Perform the final assembly
  void finalAssembly();

  /// Start to assemble the matrix
  virtual void beginAssembly(LSSMatrixAssemblyType assemblyType) = 0;

  /// Finish to assemble the matrix
  virtual void endAssembly(LSSMatrixAssemblyType assemblyType) = 0;

  /// Print this matrix
  virtual void printToScreen() const = 0;

  /// Print this matrix to a file
  virtual void printToFile(const char* fileName) const = 0;

  /// Set all the elements of the Vector equal to the given value
  virtual void setValue(const CFint im,
                        const CFint in,
                        const CFreal value) = 0;

  /// Set a list of values
  virtual void setValues(const CFuint m,
                         const CFint* im,
                         const CFuint n,
                         const CFint* in,
                         const CFreal* values) = 0;

  /// Add an array of values in the vector
  virtual void addValue(const CFint im,
                        const CFint in,
                        const CFreal value) = 0;

  /// Add a list of values
  virtual void addValues(const CFuint m,
                         const CFint* im,
                         const CFuint n,
                         const CFint* in,
                         const CFreal* values) = 0;

  /// Get one value
  virtual void getValue(const CFint im,
                        const CFint in,
                        CFreal& value) = 0;
  /// Get a list of values
  virtual void getValues(const CFuint m,
                         const CFint* im,
                         const CFuint n,
                         const CFint* in,
                         CFreal* values) = 0;

  /// Set a row, diagonal and off-diagonals 
  virtual void setRow(const CFuint row, CFreal diagval, CFreal offdiagval) = 0;

  /// Set the diagonal
  virtual void setDiagonal(LSSVector& diag) = 0;

  /// Add to the diagonal
  virtual void addToDiagonal(LSSVector& diag) = 0;

  /// Reset to 0 all the non-zero elements of the matrix
  virtual void resetToZeroEntries() = 0;

  /// Set a list of values
  virtual void setValues(const BlockAccumulator& acc) = 0;

  /// Add a list of values
  virtual void addValues(const BlockAccumulator& acc) = 0;

  /// Freeze the matrix structure concerning the non zero
  /// locations
  virtual void freezeNonZeroStructure() = 0;

private:

  /// Copy constructor
  LSSMatrix(const LSSMatrix& other);

  /// Overloading of the assignment operator
  const LSSMatrix& operator= (const LSSMatrix& other);

protected:
  
  /// set on GPU
  bool m_useGPU;
  
}; // end of class LSSMatrix

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LSSMatrix_hh
