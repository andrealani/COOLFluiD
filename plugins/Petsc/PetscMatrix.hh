// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_PetscMatrix_hh
#define COOLFluiD_Numerics_Petsc_PetscMatrix_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "MathTools/RealVector.hh"
#include "Framework/LSSMatrix.hh"

#include "Petsc/PetscVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class BlockAccumulator; }

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a proxy for the PETSC Mat object.
 *
 * @author Andrea Lani
 * @author Jiri Simonek
 *
 */
class PetscMatrix : public Framework::LSSMatrix {

public: // functions

  /**
   * Default constructor without arguments.
   */
  PetscMatrix();

  /**
   * Destructor.
   */
  ~PetscMatrix();

  /**
   * Create a sequential sparse matrix
   * @param m number of Rows of the matrix
   * @param n number of Columns of the matrix
   * @param nz maximum number of non-zero entries
   * @param nnz exact number of non-zero entries in the diagonal portion of the local matrix for each row
   */
  void createSeqAIJ(const CFint m,
		    const CFint n,
		    const CFint nz,
		    const CFint* nnz,
		    const char* name = CFNULL);
  
  /**
   * Create a sequential block sparse matrix
   * @param blocksize size of the block matrix
   * @param m number of Rows of the matrix
   * @param n number of Columns of the matrix
   * @param nz maximum number of non-zero entries
   * @param nnz exact number of non-zero entries in the diagonal portion of the local matrix for each row
   */
  void createSeqBAIJ(const CFuint blockSize,
		     const CFint m,
		     const CFint n,
		     const CFint nz,
		     const CFint* nnz,
		     const char* name = CFNULL);
  
  /**
   * Create a sequential Jacobian-Free matrix
   * @param m number of Rows of the matrix
   * @param n number of Columns of the matrix
   * @param M number of Rows of the matrix - same as "m" - (different only in parallel)
   * @param N number of Columns of the matrix - same as "n" - (different only in parallel)
   * @param ctx pointer to data needed for the Jacobian-Free method
   */
  void createSeqJFMat(const CFint m,
		      const CFint n,
		      const CFint M,
		      const CFint N,
		      void* ctx,
		      const char* name = CFNULL);
  
  /**
   * Create a parallel sparse matrix
   * @param m localNumber of Rows of the matrix
   * @param n local Number of Columns of the matrix
   * @param M Global Number of Rows of the matrix
   * @param N Global Number of Columns of the matrix
   * @param dnz maximum number of non-zero entries in the diagonal portion of the local matrix
   * @param dnnz exact number of non-zero entries in the diagonal portion of the local matrix for each row
   * @param onz maximum number of non-zero entries in the 'out of diagonal' portion of the local matrix
   * @param onnz exact number of non-zero entries in the 'out-of-diagonal' portion of the local matrix for each row
   */
#ifdef CF_HAVE_MPI
  void createParAIJ(MPI_Comm comm,
		    const CFint m,
		    const CFint n,
		    const CFint M,
		    const CFint N,
		    const CFint dnz,
		    const CFint* dnnz,
		    const CFint onz,
		    const CFint* onnz,
		    const char* name = CFNULL);
#endif
  
  /**
   * Create a parallel block sparse matrix
   * @param blocksize size of the block matrices
   * @param m local Number of Rows of the matrix
   * @param n local Number of Columns of the matrix
   * @param M Global Number of Rows of the matrix
   * @param N Global Number of Columns of the matrix
   * @param dnz maximum number of non-zero entries in the diagonal portion of the local matrix
   * @param dnnz exact number of non-zero entries in the diagonal portion of the local matrix for each row
   * @param onz maximum number of non-zero entries in the 'out of diagonal' portion of the local matrix
   * @param onnz exact number of non-zero entries in the 'out-of-diagonal' portion of the local matrix for each row
   */
#ifdef CF_HAVE_MPI
  void createParBAIJ(MPI_Comm comm,
		     const CFuint blockSize,
		     const CFint m,
		     const CFint n,
		     const CFint M,
		     const CFint N,
		     const CFint dnz,
		     const CFint* dnnz,
		     const CFint onz,
		     const CFint* onnz,
		     const char* name = CFNULL);
#endif
  
  /**
	* Create a parallel Jacobian-Free matrix
	* @param comm MPI communicator
	* @param m local number of Rows of the matrix
	* @param n local number of Columns of the matrix
	* @param M global number of Rows of the matrix
	* @param N global number of Columns of the matrix
	* @param ctx pointer to data needed for the Jacobian-Free method
	*/
#ifdef CF_HAVE_MPI
  void createParJFMat(MPI_Comm comm,
		      const CFint m,
		      const CFint n,
		      const CFint M,
		      const CFint N,
		      void* ctx,
		      const char* name = CFNULL);
#endif
  
  /**
   * Start to assemble the matrix
   */
  void beginAssembly(LSSMatrix::LSSMatrixAssemblyType assemblyType)
  {
    MatAssemblyType matAssType = (assemblyType == FLUSH_ASSEMBLY) ?
      MAT_FLUSH_ASSEMBLY : MAT_FINAL_ASSEMBLY;
    CF_CHKERRCONTINUE(MatAssemblyBegin(m_mat, matAssType));

    // some sanity checks 
    /*MatInfo info;
    CF_CHKERRCONTINUE(MatGetInfo(m_mat,MAT_LOCAL,&info));
    std::cout << "### PetscMatrix::beginAssembly() => assembly report ###" << std::endl;
    std::cout << "### info.mallocs      = " << info.mallocs << std::endl; 
    std::cout << "### info.nz_allocated = " << info.nz_allocated << std::endl;  
    std::cout << "### info.nz_unneeded  = " << info.nz_unneeded << std::endl;    
    std::cout << "### info.nz_used      = " << info.nz_used << std::endl;
    std::cout << "### info.assemblies   = " << info.assemblies << std::endl;
    */
  }

  /**
   * Finish to assemble the matrix
   */
  void endAssembly(LSSMatrix::LSSMatrixAssemblyType assemblyType)
  {
    MatAssemblyType matAssType = (assemblyType == FLUSH_ASSEMBLY) ?
      MAT_FLUSH_ASSEMBLY : MAT_FINAL_ASSEMBLY;
    CF_CHKERRCONTINUE(MatAssemblyEnd(m_mat, matAssType));
  }

  /**
   * Print this matrix
   */
  void printToScreen() const;

  /**
   * Print this matrix to a file
   */
  void printToFile(const char* fileName) const;

  /**
   * Set all the elements of the Vector equal to the given value
   */
  void setValue(const CFint im,
                const CFint in,
                const CFreal value)
  {
    CF_CHKERRCONTINUE(MatSetValues(m_mat, 1, &im, 1, &in, &value,
        INSERT_VALUES));
  }

  /**
   * Set a list of values
   */
  void setValues(const CFuint m,
                 const CFint* im,
                 const CFuint n,
                 const CFint* in,
                 const CFreal* values)
  {
    CF_CHKERRCONTINUE(MatSetValues(m_mat, m, im, n, in, values,
                               INSERT_VALUES));
  }

  /**
   * Add an array of values in the vector
   */
  void addValue(const CFint im,
                const CFint in,
                const CFreal value)
  {
    CF_CHKERRCONTINUE(MatSetValues(m_mat, 1, &im, 1, &in, &value,
                               ADD_VALUES));
  }

  /**
   * Add a list of values
   */
  void addValues(const CFuint m,
                 const CFint* im,
                 const CFuint n,
     const CFint* in,
     const CFreal* values)
  {
    CF_CHKERRCONTINUE(MatSetValues(m_mat, m, im, n, in, values,
        ADD_VALUES));
  }

  /**
   * Get a list of values
   */
  void getValue(const CFint im,
                const CFint in,
                CFreal& value)
  {
    CF_CHKERRCONTINUE(MatGetValues(m_mat, 1, &im, 1, &in, &value));
  }

  /**
   * Get a list of values
   */
  void getValues(const CFuint m,
		 const CFint* im,
		 const CFuint n,
		 const CFint* in,
		 CFreal* values)
  {
    CF_CHKERRCONTINUE(MatGetValues(m_mat, m, im, n, in, values));
  }


  /**
   * Set a row, diagonal and off-diagonals separate values 
   */
  virtual void setRow(CFuint row, CFreal diagval, CFreal offdiagval){
    CFint ncols,*cols;
    CF_CHKERRCONTINUE(MatGetRow(m_mat,row,&ncols,(const PetscInt**)&cols,0));
    PetscScalar *val=new PetscScalar[ncols],dv=(PetscScalar)diagval,odv=(PetscScalar)offdiagval;
    for(CFint i=0; i<ncols; i++) val[i] = ((CFint)row==cols[i]) ? dv : odv;
    CF_CHKERRCONTINUE(MatSetValues(m_mat,1,(const PetscInt*)&row,ncols,cols,val,INSERT_VALUES));
    finalAssembly();
    CF_CHKERRCONTINUE(MatRestoreRow(m_mat,row,&ncols,(const PetscInt**)&cols,0));
    delete(val);
  }
  
  /**
   * Set the diagonal
   */
  void setDiagonal(Framework::LSSVector& diag)
  {
    CF_CHKERRCONTINUE(MatDiagonalSet(m_mat,
      dynamic_cast<PetscVector&>(diag).getVec(),
      INSERT_VALUES));
  }

  /**
   * Add to the diagonal
   */
  void addToDiagonal(Framework::LSSVector& diag)
  {
    CF_CHKERRCONTINUE(MatDiagonalSet(m_mat,
          dynamic_cast<PetscVector&>(diag).getVec(),
          ADD_VALUES));
  }

  /**
   * Reset to 0 all the non-zero elements of the matrix
   */
  void resetToZeroEntries()
  {
    if (!_isMatShell) {
      CF_CHKERRCONTINUE(MatZeroEntries(m_mat));
    }
  }

  /**
   * Set a list of values
   */
  void setValues(const Framework::BlockAccumulator& acc);

  /**
   * Add a list of values
   */
  void addValues(const Framework::BlockAccumulator& acc);

  /**
   * Freeze the matrix structure concerning the non zero
   * locations
   */
  void freezeNonZeroStructure()
  {
    MatSetOption(m_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
   // MatSetOption(m_mat, MAT_NEW_NONZERO_LOCATIONS_ERR, PETSC_FALSE);
  }

  /**
   * Gets the matrix
   */
  Mat getMat() const
  {
    return m_mat;
  }
	
  /**
   * Set Jacobian-Free function which provides Jacobian finite difference approach
   * @param JFFunction pointer to this function
   */
  void setJFFunction(void (*JFFunction)(void));
  
  /**
   * Get Jacobian-Free context
   * @param ctx JF method context
   */
  void getJFContext(void **ctx);
  
  /**
   * Set the AIJ flag
   * @param isAIJ flag
   */
  void setAIJ(bool isAIJ);
  
private: // data

  /// matrix
  Mat m_mat;

  /// flag to tell if the matrix is a shell matrix
  bool _isMatShell;
  
  /// flag to tell if the matrix is a AIJ
  bool _isAIJ;
  
}; // end of class PetscMatrix

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_PetscMatrix_hh
