// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_TrilinosMatrix_hh
#define COOLFluiD_Numerics_Trilinos_TrilinosMatrix_hh

//////////////////////////////////////////////////////////////////////////////




#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include "Trilinos.hh"
#include "TrilinosVector.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "MathTools/RealVector.hh"


#define linefile printf("Line %u in file %s.\n",__LINE__,__FILE__);fflush(stdout);

//#define ADD_ACCU_ROW_BY_ROW  // TODO --> must sure that accu.getWorkspace() is correctly allocated

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a proxy for an Epetra matrix.
 *
 * PROCEDURE:
 * -# construction
 * -# set the Epetra_Map
 * -# creation: size arguments are redundant because already available in the Epetra_Map
 *              (This cannot be avoided, because deriving from LSSMatrix ...)
 * -# building the matrix: set/add-routines, with rules:
 *    - before calling begin/endAssembly(), anything can be set/added to any location in the matrix
 *    - after calling begin/endAssembly(), only set/add to existing (= nonzero) entries will have an effect.
 *      So, setting/adding to a zero location simply has NO effect at all.
 *      If NDEBUG is not defined, warnings will be issued when this happens.
 * -# beginAssembly()
 * -# endAssembly(): construction of the nonzero structure
 * -# changing existing entries (see also step 4.)
 *
 * RESTART FROM SCRATCH: destroy()
 * - destroying the pointer to the map, destroying the underlying matrix
 * - afterwards, the above procedure needs to be restarted from step 2 (setting the Epetra_Map)
 *
 * RESTART REUSING THE EXISTING NONZERO STRUCTURE:
 * -# freezeNonzeroStructure() or begin/endAssembly()
 * -# resetToZeroEntries()
 * - the above procedure can be restarted from 7 (changing existing entries)
 *
 * @author Tim Boonen
 *
 */
class TrilinosMatrix : public Framework::LSSMatrix {
public:

  /**
   * Default constructor without arguments.
   */
  TrilinosMatrix();

  /**
   * Destructor.
   */
  ~TrilinosMatrix();

  /**
   *  Deallocates/deletes all internal data: deletion of the pointer
   *  to the map (the map itself is not deleted), and deletion of all other
   *  matrix data (including the nonzero structure) except the name.
   */
  void destroy();

  /**
   * Sets the Epetra map for the unknowns (not for the states !!)
   */
  void setEpetraMap(const Epetra_Map *map)
  {
    _map = map;
  }

  /**
   *  Returns the Epetra map for the unknowns (not for the states !!)
   */
  const Epetra_Map* getEpetraMap()
  {
    return _map;
  }

  /**
   * Create a sequential sparse matrix
   * @param     m       number of rows (individual, NOT state)
   * @param     n       number of columns (individual, NOT state)
   * @param     nz      number of nonzero entries per row (equal for all rows)
   * @param     nnz     array containing the number of nonzero entries in the rows
   *                    (possibly different for each row) or NULL
   */
  void createSeqAIJ(const CFint m,
        const CFint n,
        const CFint nz,
        const int* nnz,
        const char* name = "");

  /**
   * Create a sequential block sparse matrix
   * @param     m       number of individual rows (NOT block rows !!!)
   * @param     n       number of individual columns (NOT block columns !!!)
   * @param     nz      number of nonzero block entries per BLOCK row (equal for all block rows)
   * @param     nnz     array containing the number of nonzero BLOCKS in the various BLOCK rows
   *                    (possibly different for each block row) or PETSC_NULL
   * @remarks   The underlying Epetra matrix is NOT a block matrix at this moment.
   */
  void createSeqBAIJ(const CFuint blockSize,
         const CFint m,
         const CFint n,
         const CFint nz,
         const int* nnz,
         const char* name = "");

  /**
   * Create a parallel sparse matrix
   * @param     m       local number of individual rows
   * @param     n       local number of individual columns
   * @param     M       global number of individual rows
   * @param     N       global number of individual columns
   */
#ifdef CF_HAVE_MPI
  void createParAIJ(MPI_Comm comm,
        const CFint m,
        const CFint n,
        const CFint M,
        const CFint N,
        const CFint dnz,
        const int* dnnz,
        const CFint onz,
        const int* onnz,
        const char* name = "");
#endif

  /**
   * Create a parallel block sparse matrix
   * @param     m       local number of individual rows (NOT row states !!!)
   * @param     n       local number of individual columns (NOT coolumn states !!!)
   * @param     M       global number of individual rows (NOT row states !!!)
   * @param     N       global number of individual columns (NOT columns states !!!)
   * @param     nz      number of nonzero BLOCK entries per BLOCK row (equal for all block rows)
   * @param     nnz     array containing the number of nonzero BLOCKS in the various BLOCK rows
   *                    (possibly different for each block row) or PETSC_NULL
   * @remarks   The underlying Epetra matrix is NOT a block matrix at this moment.
   */
#ifdef CF_HAVE_MPI
  void createParBAIJ(MPI_Comm comm,
         const CFuint blockSize,
         const CFint m,
         const CFint n,
         const CFint M,
         const CFint N,
         const CFint dnz,
         const int* dnnz,
         const CFint onz,
         const int* onnz,
         const char* name = "");
#endif

  /**
   * Start to assemble the matrix
   */
  void beginAssembly(COOLFluiD::Framework::LSSMatrix::LSSMatrixAssemblyType type=COOLFluiD::Framework::LSSMatrix::FINAL_ASSEMBLY)
  {
    cf_assert(_toBeDestroyed == true);
    // nothing to do for trilinos
  }

  /**
   * Finish to assemble the matrix
   */
  void endAssembly(COOLFluiD::Framework::LSSMatrix::LSSMatrixAssemblyType type=COOLFluiD::Framework::LSSMatrix::FINAL_ASSEMBLY)
  {
    cf_assert(_toBeDestroyed == true);
    if (!isAssembled()) {
      _mat->FillComplete();
      _mat->OptimizeStorage();
      _isAssembled = true;
    }
  }

  /**
   * Print this matrix to screen
   */
  void printToScreen() const
  {
    cf_assert(isAssembled()); // otherwise, the view is not always clear (duplicate entries ...)
    CFout << "TRILINOSMATRIX " << _name << "\n";
    CFout << "************** " << "\n";
    CFout << *_mat << "\n";
  }

  /**
   * Print this matrix to a file (NOT IMPLEMENTED)
   */
  void printToFile(const char* fileName) const;

  /**
   *  Builds the finite difference discretization of
   *  1D poisson equation. (A(i,i)=2, A(i,i-1)=-1=A(i,i+1))
   *  @remarks  it is not needed to call begin/endAssembly() afterwards
   */
  void set1DPoisson(const Epetra_Map *map);

  /**
   * Set an element of the Matrix equal to the given value
   * @param     row     global row (NOT block row !!!)
   * @param     col     global column (NOT block column !!!)
   * @remarks   if row is not owned by this processor, nothing will be done !!!
   */
  void setValue(CFint row,
                CFint col,
                CFreal valueToSet)
  {
    // checking
    cf_assert(_toBeDestroyed == true);

    if (_mat->MyGRID(row)) {
      if (isAssembled()) {
        cf_assert(_mat->IndicesAreLocal());
        int lrid = _mat->LRID(row);
        int lcid = _mat->LCID(col);
        int ierr = _mat->ReplaceMyValues(lrid, 1, &valueToSet, &lcid);
        cf_assert(ierr==0);  // if not, then the entry to be replaced did not exist
                          //         --> replacement did not succeed
      }
      else {
         int ierrElse = _mat->ReplaceGlobalValues(row, 1, &valueToSet, &col);            // ierrElse>0: entry not present yet
         if (ierrElse>0) ierrElse = _mat->InsertGlobalValues(row, 1, &valueToSet, &col); // ierrElse>0: row allocation not sufficient
         cf_assert(ierrElse==0); // if ierrElse>0, preallocated memory was not sufficient.
      }
    }
  }

  /**
   * Set a list of values using global indices of the individual unknowns
   * (not of the states !!!)
   * @TODO make more efficient
   */
  void setValues(CFuint nbRows,
                 const CFint* rows,
                 CFuint nbCols,
                 const CFint* cols,
                 const CFreal* values)
  {
    for (CFuint r=0; r<nbRows; r++) {
      for (CFuint c=0; c<nbCols; c++) {
        setValue(rows[r], cols[c], values[r*nbCols + c]);
      }
    }
  }

  /**
   * Add an array of values in the vector
   * @param     row     global row (NOT block row !!!)
   * @param     col     global column (NOT block column !!!)
   * @remarks   if row is not owned by this processor, nothing will be done !!!
   */
  void addValue(CFint row,
                CFint col,
                CFreal valueToAdd)
  {
    // checking
    cf_assert(_toBeDestroyed == true);

    if (_mat->MyGRID(row)) {
      if (isAssembled()) {
        cf_assert(_mat->IndicesAreLocal());
        int lrid = _mat->LRID(row);
        int lcid = _mat->LCID(col);
        int ierr = _mat->SumIntoMyValues(lrid, 1, &valueToAdd, &lcid);
        cf_assert(ierr == 0);      // fails if replacing did not succeed
                                // (replacing nonexisting entry has no effect in epetra)
      }
      else {
         int ierrElse = _mat->SumIntoGlobalValues(row, 1, &valueToAdd, &col);            // ierrElse>0: entry not present yet
         if (ierrElse>0) ierrElse = _mat->InsertGlobalValues(row, 1, &valueToAdd, &col); // ierrElse>0: row allocation not sufficient
         cf_assert(ierrElse==0);
      }
    }
  }


  /**
   * Add a list of values using global indices of the individual unknowns
   * (not of the states !!!)
   * @TODO make more efficient
   */
  void addValues(CFuint nbRows,
                 const CFint* rows,
                 CFuint nbCols,
                 const CFint* cols,
                 const CFreal* values)
  {
    cf_assert(_toBeDestroyed == true);
    for (CFuint r=0; r<nbRows; r++) {
      for (CFuint c=0; c<nbCols; c++) {
        addValue(rows[r], cols[c], values[r*nbCols + c]);
      }
    }
  }

  /**
   * Returns an entry
   * @param     row     global row index (of individual unknown, not of state !!!)
   * @param     col     global column index (of individual unknown, not of state !!!)
   * @pre       row must be locally owned
   */
  double getValue(int row, int col)
  {
    cf_assert(isAssembled());
    //cf_assert(_mat->MyGRID(row));

    if (!_mat->MyGRID(row)) {
      CFerr << "WARNING: TrilinosMatrix::getValue(row,col) for nonlocal row --> simply 0 returned" << "\n";
      return 0.0;
    }

    int    nb, *index;
    double *val;
    int    c;
    if (_mat->IndicesAreLocal()) {
      int rowLID = _mat->LRID(row);
      c = _mat->LCID(col);
      _mat->ExtractMyRowView(rowLID, nb, val, index);
    }
    else {
      c = col;
      _mat->ExtractGlobalRowView(row, nb, val, index);
    }
    for (int i=0; i<nb; i++) {
      if (index[i] == c) {
        return val[i];
      }
    }
    return 0.0;
  }

  /**
   * Get a list of values according to global indices of individual rows and columns
   * (the values are ordered row by row)
   * @pre       all rows must be locally owned
   */
  void getValues(CFuint nbRows,
                 const CFint* rows,
                 CFuint nbCols,
                 const CFint* cols,
                 CFreal* values)
  {
    cf_assert(isAssembled());
    for (CFuint r=0; r<nbRows; r++) {
      for (CFuint c=0; c<nbCols; c++) {
        values[r*nbCols + c] = getValue(rows[r], cols[c]);
      }
    }
  }

  /**
   * Get one value
   */
  void getValue(const CFint im,
                const CFint in,
                CFreal& value)
  {
    value = getValue(im, in);
  }

  /**
   * Set a row, diagonal and off-diagonals separate values 
   */
  virtual void setRow(const CFuint row, CFreal diagval, CFreal offdiagval){
    int numentries,*indices;
    double *values,dv=(double)diagval,odv=(double)offdiagval;
    int rowLID = _mat->LRID(row);
    int err= _mat->ExtractMyRowView(rowLID,numentries,values,indices);
    for(CFuint i=0; i<numentries; i++) values[i]= rowLID==indices[i] ? dv : odv;
  }

  /**
   * Set the diagonal
   */
  void setDiagonal(COOLFluiD::Framework::LSSVector& diag)
  {
    CFerr << "TrilinosMatrix::setDiagonal(TrilinosVector&) --> not implemented" << "\n";
    exit(1);
  }

  /**
   * Add to the diagonal
   */
  void addToDiagonal(COOLFluiD::Framework::LSSVector& diag)
  {
    CFerr << "TrilinosMatrix::setDiagonal(TrilinosVector&) --> not implemented" << "\n";
    exit(1);
  }

  /**
   * Reset to 0 all the non-zero elements of the matrix
   */
  void resetToZeroEntries()
  {
    //cf_assert(isAssembled()); // REMARK: this requirement is apparently stricter than the way coolfluid uses the LSS-interface ...
    cf_assert(_toBeDestroyed == true);
    _mat->PutScalar(0.0);
  }

  /**
   * Set a list of values
   * TODO make more efficient
   */
void setValues(const Framework::BlockAccumulator& accu)
{
  cf_assert(_toBeDestroyed == true);

#ifdef ADD_ACCU_ROW_BY_ROW

  // fetching information from the accumulator
  int nbBlockRows = accu.getM();
  int nbBlockCols = accu.getN();
  int blockSize = accu.getNB();
  int nbRows = nbBlockRows * blockSize;
  int nbCols = nbBlockCols * blockSize;

  cf_assert((int)accu.getWorkspace().size() >= (int)blockSize*(nbBlockRows+nbBlockCols));  // otherwise, the workspace is not correctly allocated !!

  // calculation of the global indices
  int *rid = const_cast<int*>(&(accu.getWorkspace()[0]));       // TODO: avoid const-casting
  int *cid = const_cast<int*>(&(accu.getWorkspace()[nbRows]));  // TODO: avoid const-casting
  for (int R=0; R<nbBlockRows; R++) {
    int offset = accu.im[R] * blockSize;
    for (int r=0; r<blockSize; r++) {
      rid[blockSize*R+r] = offset + r;
    }
  }
  for (int C=0; C<nbBlockCols; C++) {
    int offset = accu.getIN()[C] * blockSize;
    for (int c=0; c<blockSize; c++) {
      cid[blockSize*C+c] = offset + c;
    }
  }

  // setting
  double *val = const_cast<double*>(&accu.values[0]);
  if (_mat->IndicesAreLocal()) {
    // global to local ID translation
    for (int r=0; r<nbRows; r++) rid[r] = _mat->LRID(rid[r]);
    for (int c=0; c<nbCols; c++) cid[c] = _mat->LCID(cid[c]);
    for (int r=0; r<nbRows; r++) {
      int ierr = _mat->ReplaceMyValues(rid[r], nbCols, &val[nbCols*r], cid);
      if (ierr!=0) CFout << "            setValues(accu) --> ReplaceMyValues ierr = " << ierr << "\n";
      cf_assert(ierr == 0);
    }
  }
  else {
    /* NEWER VERSION: appears to be slower ...
       bool found;
       for (int r=0; r<nbRows; r++) {

       // declarations
       int    nb, *found_col;
       double *found_val;

       // fetching the existing entries
       _mat->ExtractGlobalRowView(rid[r], nb, found_val, found_col);

       // setting
       for (int c=0; c<nbCols; c++) {
       found = false;
       for (int i=0; i<nb; i++) {
       if (found_col[i] == cid[c]) {
       found = true;
       found_val[i] = val[nbCols*r+c];
       break;
       }
       }
       if (!found) {
       int ierr = _mat->InsertGlobalValues(rid[r], 1, &val[nbCols*r+c], &cid[c]);
       if (ierr!=0) CFout << "     setValues(accu) --> InsertGlobalValues ierr = " << ierr << "\n";
       }
       }
       }*/

    // OLD VERSION: uses unneeded Insert-routines, hence requiring >preallocated memory per row
    for (int r=0; r<nbRows; r++) {
      //_mat->InsertGlobalValues(rid[r], nbCols, &val[nbCols*r], cid);
      _mat->InsertGlobalValues(rid[r], nbCols, &val[nbCols*r], cid);
      //if (ierr!=0) CFout << "            setValues(accu) --> InsertGlobalValues ierr = " << ierr << "\n";
      // @remarks: resulting ierr>0 = preallocated nnz not sufficient !!!
    }
  }


#else // ADD_ACCU_ROW_BY_ROW

  // fetching information from the accumulator
  int nbBlockRows = accu.getM();
  int nbBlockCols = accu.getN();
  int blockSize = accu.getNB();
  int nbRows = nbBlockRows * blockSize;
  int nbCols = nbBlockCols * blockSize;

  // translation of the indices
  int *row = new int[nbBlockRows * blockSize];
  int *col = new int[nbBlockCols * blockSize];
  for (int R=0; R<nbBlockRows; R++) {
    int offset = accu.getIM()[R] * blockSize;
    for (int r=0; r<blockSize; r++) {
      row[blockSize*R+r] = offset + r;
    }
  }
  //CFout << "global row indices : ";
  //for (int i=0; i<nbRows; i++) CFout << row[i] << " ";
  //CFout << "\n";
  for (int C=0; C<nbBlockCols; C++) {
    int offset = accu.getIN()[C] * blockSize;
    for (int c=0; c<blockSize; c++) {
      col[blockSize*C+c] = offset + c;
    }
  }
  //CFout << "global col indices : ";
  //for (int i=0; i<nbCols; i++) CFout << col[i] << " ";
  //CFout << "\n";

  // setting: TODO --> row by row = more efficient
  for (int r=0; r<nbRows; r++) {
    for (int c=0; c<nbCols; c++) {
      setValue(row[r], col[c], ((Framework::BlockAccumulator*)(&accu))->getPtr()[nbCols*r+c]);
    }
  }
  
  // deallocation
  delete[] row;
  delete[] col;

#endif // ADD_ACCU_ROW_BY_ROW

}

  /**
   * Add a list of values
   * TODO make more efficient
   */
  void addValues(const Framework::BlockAccumulator& accu)
  {

    cf_assert(_toBeDestroyed == true);
#define ADD_ACCU_ROW_BY_ROW
#ifdef ADD_ACCU_ROW_BY_ROW

    // fetching information from the accumulator
    int nbBlockRows = accu.getM();
    int nbBlockCols = accu.getN();
    int blockSize = accu.getNB();
    int nbRows = nbBlockRows * blockSize;
    int nbCols = nbBlockCols * blockSize;

    cf_assert((int)accu.getWorkspace().size() >= (int)blockSize*(nbBlockRows+nbBlockCols));  // otherwise, the workspace is not correctly allocated !!

    // calculation of the global indices
    int *rid = const_cast<int*>(&(accu.getWorkspace()[0]));       // TODO: avoid const-casting
    int *cid = const_cast<int*>(&(accu.getWorkspace()[nbRows]));  // TODO: avoid const-casting
    for (int R=0; R<nbBlockRows; R++) {
      int offset = accu.getIM()[R] * blockSize;
      for (int r=0; r<blockSize; r++) {
        rid[blockSize*R+r] = offset + r;
      }
    }
    for (int C=0; C<nbBlockCols; C++) {
      int offset = accu.getIN()[C] * blockSize;
      for (int c=0; c<blockSize; c++) {
        cid[blockSize*C+c] = offset + c;
      }
    }

    // setting
    //    double *val = const_cast<double*>(&accu.values[0]);
    //    double *val = accu.values.ptr();
    //    RealVector rvval(accu.values); /// @warning heavy, do direct castingto tthe raw pointer of accu.values
    double *val = ((Framework::BlockAccumulator*)(&accu))->getPtr();

    if (_mat->IndicesAreLocal()) {
      cf_assert(isAssembled());
      // global to local ID translation
      for (int r=0; r<nbRows; r++) {
        if (_mat->MyGRID(rid[r])) { rid[r] = _mat->LRID(rid[r]); }
        else { rid[r]=-1; }
      }
      for (int c=0; c<nbCols; c++) cid[c] = _mat->LCID(cid[c]);
      for (int r=0; r<nbRows; r++) {
        if (rid[r]!=-1) {
          int ierr = _mat->SumIntoMyValues(rid[r], nbCols, &val[nbCols*r], cid);
          if (ierr!=0) CFout << "            addValues(accu) --> SumIntoMyValues ierr = " << ierr << "\n";
          cf_assert(ierr == 0);
        }
      }
    }
    else {
      /* NEWER VERSION: appears to be slower ...
      bool found;
      for (int r=0; r<nbRows; r++) {

        // declarations
        int    nb, *found_col;
        double *found_val;

        // fetching the existing entries
        _mat->ExtractGlobalRowView(rid[r], nb, found_val, found_col);

        // setting
        for (int c=0; c<nbCols; c++) {
          found = false;
          for (int i=0; i<nb; i++) {
            if (found_col[i] == cid[c]) {
              found = true;
              found_val[i] += val[nbCols*r+c];
              break;
            }
          }
          if (!found) {
            int ierr = _mat->InsertGlobalValues(rid[r], 1, &val[nbCols*r+c], &cid[c]);
            if (ierr!=0) CFout << "     addValues(accu) --> InsertGlobalValues ierr = " << ierr << "\n";
          }
        }

      }*/

      // OLD VERSION: uses unneeded Insert-routines, hence requiring >preallocated memory per row
      cf_assert(!isAssembled());
      for (int r=0; r<nbRows; r++) {
        //_mat->InsertGlobalValues(rid[r], nbCols, &val[nbCols*r], cid);
        _mat->InsertGlobalValues(rid[r], nbCols, &val[nbCols*r], cid);
        //if (ierr!=0) CFout << "            addValues(accu) --> InsertGlobalValues ierr = " << ierr << "\n";
        // @remarks: resulting ierr>0 = preallocated nnz not sufficient !!!
      }
    }

#else // ADD_ACCU_ROW_BY_ROW

    // fetching information from the accumulator
    int nbBlockRows = accu.getM();
    int nbBlockCols = accu.getN();
    int blockSize = accu.getNB();
    int nbRows = nbBlockRows * blockSize;
    int nbCols = nbBlockCols * blockSize;

    // translation of the indices
    int *row = new int[nbBlockRows * blockSize];
    int *col = new int[nbBlockCols * blockSize];
    for (int R=0; R<nbBlockRows; R++) {
      int offset = accu.getIM()[R] * blockSize;
      for (int r=0; r<blockSize; r++) {
        row[blockSize*R+r] = offset + r;
      }
    }
    //CFout << "global row indices : ";
    //for (int i=0; i<nbRows; i++) CFout << row[i] << " ";
    //CFout << "\n";
    for (int C=0; C<nbBlockCols; C++) {
      int offset = accu.getIN()[C] * blockSize;
      for (int c=0; c<blockSize; c++) {
        col[blockSize*C+c] = offset + c;
      }
    }
    //CFout << "global col indices : ";
    //for (int i=0; i<nbCols; i++) CFout << col[i] << " ";
    //CFout << "\n";


    // setting: TODO --> row by row = more efficient
    for (int r=0; r<nbRows; r++) {
      for (int c=0; c<nbCols; c++) {
        addValue(row[r], col[c], accu.values[nbCols*r+c]);
      }
    }

    // deallocation
    delete[] row;
    delete[] col;

#endif // ADD_ACCU_ROW_BY_ROW
#undef ADD_ACCU_ROW_BY_ROW
  }

  /**
   * Freeze the matrix structure concerning the non zero
   * locations
   */
  void freezeNonZeroStructure()
  {
    cf_assert(_toBeDestroyed);
    beginAssembly();
    endAssembly();
  }

  /**
   * Returns the underlying Epetra matrix
   */
  Epetra_CrsMatrix* getMat() const
  {
    return _mat;
  }

  /**
   *  Returns wether the matrix is assembled.
   */
  bool isAssembled() const {return _isAssembled;}


 private:

  /**
   * Copy constructor.
   * Deep copy for the Epetra_Matrix, but not for the Epetra_Map and the name.
   */
  TrilinosMatrix(const TrilinosMatrix& other);

  /**
   * Assignment operator
   * Deep copy for the Epetra_Matrix, but not for the Epetra_Map and the name.
   */
  const TrilinosMatrix& operator= (const TrilinosMatrix& other);

  /// matrix
  Epetra_CrsMatrix* _mat;

  /// Epetra map for the individual unknowns (NOT for the states !!!)
  //  deletion is NOT the task of this class
  const Epetra_Map* _map;

  /// true after creation, false after destroy()
  bool _toBeDestroyed;

  /// true after endAssembly()
  bool _isAssembled;    // if true, the nonzero structure cannot be changed any more

  /// name of the matrix
  string _name;

}; // end of class TrilinosMatrix

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_TrilinosMatrix_hh

