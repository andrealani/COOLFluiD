// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BlockAccumulator_hh
#define COOLFluiD_Framework_BlockAccumulator_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "LSSIdxMapping.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an accumulator of values to be put in a block matrix
/// @author Andrea Lani
class Framework_API BlockAccumulator {

public:

  /// Constructor
  BlockAccumulator(const CFuint nbRows,
		   const CFuint nbCols,
		   const CFuint subBlockSize,
		   const LSSIdxMapping& localToGlobal,
		   CFreal* ptr = CFNULL);
  
  /// Destructor.
  ~BlockAccumulator();
 
  /// Reset the array and all sizes
  void resetB(const CFuint nbRows,
              const CFuint nbCols,
              const CFuint subBlockSize,
              CFreal* array);
 
  /// Set a single value in all the entries
  void setValue(const CFreal value)
  {
    values = value;
  }

  /// Set a single value
  void setValue(const CFuint i,
		const CFuint j,
		const CFuint ib,
		const CFuint jb,
		const CFreal value)
  {
    const CFuint idx = getIndex(i,j,ib,jb);
    cf_assert(idx < m*n*nb*nb);
    values[idx] = value;
  }
  
  /// Set a list of values
  void setValues(const CFuint i,
		 const CFuint j,
		 const CFuint jb,
		 const CFreal* inValues)
  {
    for (CFuint ib = 0; ib < nb; ++ib) {
      const CFuint idx = getIndex(i,j,ib,jb);
      cf_assert(idx < m*n*nb*nb);
      values[idx] = inValues[ib];
    }
  }
  
  /// Add a single value in all the entries
  void addValue(const CFreal value)
  {
    values += value;
  }
  
  /// Add a single value in the vector
  void addValue(const CFuint i,
		const CFuint j,
		const CFuint ib,
		const CFuint jb,
		const CFreal value)
  {
    const CFuint idx = getIndex(i,j,ib,jb);
    cf_assert(idx < m*n*nb*nb);
    values[idx] += value;
  }
  
  /// Add an array of values in the vector
  void addValues(const CFuint i,
		 const CFuint j,
		 const CFuint jb,
		 const CFreal* inValues)
  {
    for (CFuint ib = 0; ib < nb; ++ib) {
      const CFuint idx = getIndex(i,j,ib,jb);
      cf_assert(idx < m*n*nb*nb);
      values[idx] += inValues[ib];
    }
  }
  
  /// Sets a RealMatrix of values in the vector
  void setValuesM(const CFuint i,
		  const CFuint j,
		  const RealMatrix& matrix)
  {
    cf_assert(matrix.nbRows() == nb);
    cf_assert(matrix.nbCols() == nb);
    for (CFuint ib = 0; ib < nb; ++ib) {
      for (CFuint jb = 0; jb < nb; ++jb) {
	const CFuint idx = getIndex(i,j,ib,jb);
	cf_assert(idx < m*n*nb*nb);
	values[idx] = matrix(ib,jb);
      }
    }
  }
  
  /// Copy the entries of a matrix in the vector
  void setValuesM(const RealMatrix& matrix)
  {
    cf_assert(matrix.nbRows() == m*nb);
    cf_assert(matrix.nbCols() == n*nb);
    CFuint matCol;
    CFuint matRow;
    
    for (CFuint R=0; R<m; R++) {
      for (CFuint r=0; r<nb; r++) {
	matRow = R*nb + r;
	for (CFuint C=0; C<n; C++) {
	  for (CFuint c=0; c<nb; c++) {
	    matCol = C*nb + c;
	    values[getIndex(R,C,r,c)] = matrix(matRow,matCol);
	  }
	}
      }
    }
  }
  
  /// Adds the entries of a matrix in the vector
  void addValuesM(const RealMatrix& matrix)
  {
    cf_assert(matrix.nbRows() == m*nb);
    cf_assert(matrix.nbCols() == n*nb);
    CFuint matCol;
    CFuint matRow;
    
    for (CFuint R=0; R<m; R++) {
      for (CFuint r=0; r<nb; r++) {
	matRow = R*nb + r;
	for (CFuint C=0; C<n; C++) {
	  for (CFuint c=0; c<nb; c++) {
	    matCol = C*nb + c;
	    values[getIndex(R,C,r,c)] += matrix(matRow,matCol);
	  }
	}
      }
    }
  }
  
  /// Adds a RealMatrix of values in the vector
  void addValuesM(const CFuint i,
		  const CFuint j,
		  const RealMatrix& matrix)
  {
    cf_assert(matrix.nbRows() == nb);
    cf_assert(matrix.nbCols() == nb);
    for (CFuint ib = 0; ib < nb; ++ib) {
      for (CFuint jb = 0; jb < nb; ++jb) {
	const CFuint idx = getIndex(i,j,ib,jb);
	cf_assert(idx < m*n*nb*nb);
	values[idx] += matrix(ib,jb);
      }
    }
  }
  
  /// Set index of the row
  void setRowIndex(const CFint iRow, const CFint idx)
  {
    cf_assert(iRow < (CFint)m);
    //    im[iRow] = _localToGlobal.getRowID(idx);
    im[iRow] = (idx >= 0) ? _localToGlobal.getRowID(idx) : idx;
  }
  
  /// Set index of the column
  void setColIndex(const CFint iCol, const CFint idx)
  {
    cf_assert(iCol < (CFint)n);
    //    in[iCol] = _localToGlobal.getColID(idx);
    in[iCol] = (idx >= 0) ? _localToGlobal.getColID(idx) : idx;
  }
  
  /// Set both the indexes of the row and column
  void setRowColIndex(const CFint i, const CFint idx)
  {
    cf_assert(i < (CFint)n);
    cf_assert(i < (CFint)m);
    // im[i] = _localToGlobal.getRowID(idx);
    // in[i] = _localToGlobal.getColID(idx);
    im[i] = (idx >= 0) ? _localToGlobal.getRowID(idx) : idx;
    in[i] = (idx >= 0) ? _localToGlobal.getColID(idx) : idx;
  }
  
  /// Reset to the specified value (0. by default)
  void reset(CFreal value = 0.0) {values = value;}
  
  void printToScreen() const;

  ///  Sets the coefficients for the calculation of the index
  ///  Use only if you know what you are doing !!!
  ///  formula: values[a*i+b*j+c*ib+d*jb] = accu(i,j,ib,jb)
  ///  for row by row storage:
  ///     a = subBlockSize^2 * nbBlockCols
  ///     b = subBlockSize
  ///     c = subBlockSize * nbBlockCols
  ///     d = 1
  void setCoefIndex(CFuint a, CFuint b, CFuint c, CFuint d)
  {
    _coefR = a;
    _coefC = b;
    _coefr = c;
    _coefc = d;
  }

  /// Get the value corresponding to the given indices
  CFreal getValue(CFuint i, CFuint j, CFuint ib, CFuint jb) const
  {
    return values[getIndex(i,j,ib,jb)];
  }

  /// Print the block accumulator
  void print() const
  {
    // loop over row index of sub-blocks
    for (CFuint i = 0; i < m; ++i) {
      // loop over column index of sub-blocks
      for (CFuint j = 0; j < n; ++j) {
	std::cout << "(" << i << "," << j << ")"<< std::endl;
	for (CFuint ib = 0; ib < nb; ++ib) {
	  for (CFuint jb = 0; jb < nb; ++jb) {
	    const CFuint idx = getIndex(i,j,ib,jb);
	    cf_assert(idx < m*n*nb*nb);
	    std::cout.precision(12);  std::cout << values[idx] << " ";
	  }
	  std::cout << std::endl;
	}
      }
    }
    std::cout << std::endl;
  }
  
  /// size of the block array
  CFuint size() const {return m_bsize;}
  
  /// number of rows
  CFuint getM() const {return m;}
  
  /// number of columns
  CFuint getN() const {return n;}
  
  /// size of the subblock 
  /// @pre number of rows = number of cols for subblocks
  CFuint getNB() const {return nb;}
  
  /// Get the array pointer 
  CFreal* getPtr() {return &values[0];}
  
  /// get the row indexes
  const std::vector<CFint>& getIM() const {return im;}
  
  /// get the column indexes
  const std::vector<CFint>& getIN() const {return in;}
  
  /// get the workspace
  const std::vector<CFint>& getWorkspace() const {return workspace;}
  
  /// get the workspace
  std::vector<CFint>& getWorkspaceNC() {return workspace;}
  
private: // helper functions

  /// Get the global (in this block) index corresponding to the
  /// given i and j
  CFuint getIndex(const CFuint i,
		  const CFuint j,
		  const CFuint ib,
		  const CFuint jb) const
  {
    //return nb*((i*nb + ib)*n + j) + jb;
    return _coefR*i + _coefC*j + _coefr*ib + _coefc*jb;
  }
 
private:

  /// idx mapping object
  const LSSIdxMapping& _localToGlobal;
  
private:
  
  /// row size
  CFuint m;
  
  /// column size
  CFuint n;
  
  /// sub block size (typically == number of equations)
  CFuint nb;
  
  /// block size
  CFuint m_bsize;
  
  /// array storing the indexes of the rows
  std::vector<CFint> im;
  
  /// array storing the indexes of the columns
  std::vector<CFint> in;
  
  /// storage for the values in this block
  RealVector values;
  
  /// coefficients used in getIndex
  CFuint _coefR;
  CFuint _coefC;
  CFuint _coefr;
  CFuint _coefc;
  
  /// added by TB: workspace (for translation of indices)
  std::vector<CFint> workspace;
  
}; // end of class BlockAccumulator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BlockAccumulator_hh
