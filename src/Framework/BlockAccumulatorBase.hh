// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BlockAccumulatorBase_hh
#define COOLFluiD_Framework_BlockAccumulatorBase_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Framework/Framework.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////
    
/// This class represents an accumulator of values to be put in a block matrix
/// @author Andrea Lani
class Framework_API BlockAccumulatorBase {

public:

  /// Constructor
  BlockAccumulatorBase(const CFuint nbRows,
		       const CFuint nbCols,
		       const CFuint subBlockSize,
		       CFreal* array);
  
  /// Destructor.
  virtual ~BlockAccumulatorBase();
  
  /// size of the block array
  CFuint size() const {return m_bsize;}
  
  /// number of rows
  CFuint getM() const {return m_m;}
  
  /// number of columns
  CFuint getN() const {return m_n;}
  
  /// size of the subblock 
  /// @pre number of rows = number of cols for subblocks
  CFuint getNB() const {return m_nb;}
  
  /// Get the array pointer 
  CFreal* getPtr() const {return m_values;}
  
  /// Set a single value in all the entries
  void setValue(const CFreal value) 
  {
    for (CFuint i = 0; i < m_bsize; ++i) m_values[i] = value;
  }
  
  /// Set a single value
  void setValue(const CFuint i, 
		const CFuint j,
		const CFuint ib, 
		const CFuint jb,
		const CFreal value)
  {
    const CFuint idx = getIndex(i,j,ib,jb);
    cf_assert(idx < m_bsize);
    m_values[idx] = value;
  }
  
  /// Set a list of values
  void setValues(const CFuint i,
		 const CFuint j,
		 const CFuint jb,
		 const CFreal* inValues)
  {
    for (CFuint ib = 0; ib < m_nb; ++ib) {
      const CFuint idx = getIndex(i,j,ib,jb);
      cf_assert(idx < m_bsize);
      m_values[idx] = inValues[ib];
    }
  }
  
  /// Add a single value in all the entries
  void addValue(const CFreal value) {for (CFuint i = 0; i < m_bsize; ++i) m_values[i] += value;}
  
  /// Add a single value in the vector
  void addValue(const CFuint i,
		const CFuint j,
		const CFuint ib,
		const CFuint jb,
		const CFreal value)
  {
    const CFuint idx = getIndex(i,j,ib,jb);
    cf_assert(idx < m_bsize);
    m_values[idx] += value;
  }
  
  /// Add an array of values in the vector
  void addValues(const CFuint i,
		 const CFuint j,
		 const CFuint jb,
		 const CFreal* inValues)
  {
    for (CFuint ib = 0; ib < m_nb; ++ib) {
      const CFuint idx = getIndex(i,j,ib,jb);
      cf_assert(idx < m_bsize);
      m_values[idx] += inValues[ib];
    }
  }

  /// Adds the entries of a matrix in the vector
  void addValuesM(const RealMatrix& matrix)
  {
    cf_assert(matrix.nbRows() == m_m*m_nb);
    cf_assert(matrix.nbCols() == m_n*m_nb);
    CFuint matCol = 0;
    CFuint matRow = 0;
    
    for (CFuint R=0; R<m_m; R++) {
      for (CFuint r=0; r<m_nb; r++) {
	matRow = R*m_nb + r;
	for (CFuint C=0; C<m_n; C++) {
	  for (CFuint c=0; c<m_nb; c++) {
	    matCol = C*m_nb + c;
	    m_values[getIndex(R,C,r,c)] += matrix(matRow,matCol);
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
    cf_assert(matrix.nbRows() == m_nb);
    cf_assert(matrix.nbCols() == m_nb);
    for (CFuint ib = 0; ib < m_nb; ++ib) {
      for (CFuint jb = 0; jb < m_nb; ++jb) {
	const CFuint idx = getIndex(i,j,ib,jb);
	cf_assert(idx < m_m*m_n*m_nb*m_nb);
	m_values[idx] += matrix(ib,jb);
      }
    }
  }
  
  /// Reset to the specified value (0. by default)
  void reset(const CFreal value = 0.0) {setValue(value);}
  
  /// Reset the array and all sizes
  void reset(const CFuint nbRows,
	     const CFuint nbCols,
	     const CFuint subBlockSize,
	     CFreal* array);
  
  ///  Sets the coefficients for the calculation of the index
  ///  Use only if you know what you are doing !!!
  ///  formula: values[a*i+b*j+c*ib+d*jb] = accu(i,j,ib,jb)
  ///  for row by row storage:
  ///     a = subBlockSize^2 * nbBlockCols
  ///     b = subBlockSize
  ///     c = subBlockSize * nbBlockCols
  ///     d = 1
  void setCoefIndex(const CFuint a, const CFuint b, 
		    const CFuint c, const CFuint d)
  {
    m_coefR = a; m_coefC = b; m_coefr = c; m_coefc = d;
  }
  
  /// Get the value corresponding to the given indices
  CFreal getValue(const CFuint i, const CFuint j, 
		  const CFuint ib, const CFuint jb) const
  {
    return m_values[getIndex(i,j,ib,jb)];
  }

  /// Print the block accumulator
  void print() const;
  
 protected: // helper functions
  
  /// Get the global (in this block) index corresponding to the
  /// given i and j
  CFuint getIndex(const CFuint i,
		  const CFuint j,
		  const CFuint ib,
		  const CFuint jb) const
  {
    //return nb*((i*nb + ib)*n + j) + jb;
    return m_coefR*i + m_coefC*j + m_coefr*ib + m_coefc*jb;
  }
  
 protected: // member data
  
  /// row size
  CFuint m_m;
  
  /// column size
  CFuint m_n;
  
  /// sub block size (typically == number of equations)
  CFuint m_nb;
  
  /// block size
  CFuint m_bsize;
  
  /// coefficients used in getIndex
  CFuint m_coefR;
  CFuint m_coefC;
  CFuint m_coefr;
  CFuint m_coefc;
  
  /// pointer to the values in this block
  CFreal* m_values;
  
}; // end of class BlockAccumulatorBase
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BlockAccumulatorBase_hh
