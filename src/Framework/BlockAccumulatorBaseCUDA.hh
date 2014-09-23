// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BlockAccumulatorBaseCUDA_hh
#define COOLFluiD_Framework_BlockAccumulatorBaseCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an accumulator of values to be put in a block matrix
/// on the device
/// @author Andrea Lani
class Framework_API BlockAccumulatorBaseCUDA {

public:

  /// Constructor
  HOST_DEVICE BlockAccumulatorBaseCUDA(const CFuint nbRows,
				       const CFuint nbCols,
				       const CFuint subBlockSize,
				       CFreal* array) : 
    m_m(nbRows), 
    m_n(nbCols),
    m_nb(subBlockSize),
    m_bsize(m_m*m_n*m_nb*m_nb) 
  { 
    // assign the array pointer
    m_values = array;
    
    // default = row by row storage
    setCoefIndex(m_nb*m_nb*m_n, m_nb, m_nb*m_n, 1);
  }
  
  /// Destructor.
  HOST_DEVICE ~BlockAccumulatorBaseCUDA() {}
  
  /// size of the block array
  HOST_DEVICE CFuint size() const {return m_bsize;}
  
  /// number of rows
  HOST_DEVICE CFuint getM() const {return m_m;}
  
  /// number of columns
  HOST_DEVICE CFuint getN() const {return m_n;}
  
  /// size of the subblock 
  /// @pre number of rows = number of cols for subblocks
  HOST_DEVICE CFuint getNB() const {return m_nb;}
  
  /// Get the array pointer 
  HOST_DEVICE CFreal* getPtr() const {return m_values;}
  
  /// Set a single value in all the entries
  HOST_DEVICE void setValue(const CFreal value) 
  {
    for (CFuint i = 0; i < m_bsize; ++i) m_values[i] = value;
  }
  
  /// Set a single value
  HOST_DEVICE void setValue(const CFuint i, 
			    const CFuint j,
			    const CFuint ib, 
			    const CFuint jb,
			    const CFreal value)
  {
    const CFuint idx = getIndex(i,j,ib,jb);
    m_values[idx] = value;
  }
  
  /// Set a list of values
  HOST_DEVICE void setValues(const CFuint i,
			     const CFuint j,
			     const CFuint jb,
			     const CFreal* inValues)
  {
    for (CFuint ib = 0; ib < m_nb; ++ib) {
      const CFuint idx = getIndex(i,j,ib,jb);
      m_values[idx] = inValues[ib];
    }
  }
  
  /// Add a single value in all the entries
  HOST_DEVICE void addValue(const CFreal value) {for (CFuint i = 0; i < m_bsize; ++i) m_values[i] += value;}
  
  /// Add a single value in the vector
  HOST_DEVICE void addValue(const CFuint i,
			    const CFuint j,
			    const CFuint ib,
			    const CFuint jb,
			    const CFreal value)
  {
    const CFuint idx = getIndex(i,j,ib,jb);
    m_values[idx] += value;
  }
  
  /// Add an array of values in the vector
  HOST_DEVICE void addValues(const CFuint i,
			     const CFuint j,
			     const CFuint jb,
			     const CFreal* inValues)
  {
    for (CFuint ib = 0; ib < m_nb; ++ib) {
      const CFuint idx = getIndex(i,j,ib,jb);
      m_values[idx] += inValues[ib];
    }
  }
  
  /// Reset to the specified value (0. by default)
  HOST_DEVICE void reset(const CFreal value = 0.0) {setValue(value);}
  
  ///  Sets the coefficients for the calculation of the index
  ///  Use only if you know what you are doing !!!
  ///  formula: values[a*i+b*j+c*ib+d*jb] = accu(i,j,ib,jb)
  ///  for row by row storage:
  ///     a = subBlockSize^2 * nbBlockCols
  ///     b = subBlockSize
  ///     c = subBlockSize * nbBlockCols
  ///     d = 1
  HOST_DEVICE void setCoefIndex(const CFuint a, const CFuint b, 
				const CFuint c, const CFuint d)
  {
    m_coefR = a; m_coefC = b; m_coefr = c; m_coefc = d;
  }
  
  /// Get the value corresponding to the given indices
  HOST_DEVICE CFreal getValue(const CFuint i, const CFuint j, 
			      const CFuint ib, const CFuint jb) const
  {
    return m_values[getIndex(i,j,ib,jb)];
  }
    
  /// Print the block accumulator
  void print() const
  {    
    // loop over row index of sub-blocks
    for (CFuint i = 0; i < m_m; ++i) {
      // loop over column index of sub-blocks
      for (CFuint j = 0; j < m_n; ++j) {
	std::cout << "(" << i << "," << j << ")"<< std::endl;
	for (CFuint ib = 0; ib < m_nb; ++ib) {
	  for (CFuint jb = 0; jb < m_nb; ++jb) {
	    const CFuint idx = getIndex(i,j,ib,jb);
	    // cf_assert(idx < m_bsize);
	    std::cout.precision(12);  std::cout << m_values[idx] << " ";
	  }
	  std::cout << std::endl;
	}
      }
    }
    std::cout << std::endl;
  }
  
 private: // helper functions
  
  /// Get the global (in this block) index corresponding to the
  /// given i and j
  HOST_DEVICE CFuint getIndex(const CFuint i,
			      const CFuint j,
			      const CFuint ib,
			      const CFuint jb) const
  {
    //return nb*((i*nb + ib)*n + j) + jb;
    return m_coefR*i + m_coefC*j + m_coefr*ib + m_coefc*jb;
  }
  
 private: // member data
  
  /// row size
  const CFuint m_m;
  
  /// column size
  const CFuint m_n;
  
  /// sub block size (typically == number of equations)
  const CFuint m_nb;
  
  /// block size
  const CFuint m_bsize;
  
  /// coefficients used in getIndex
  CFuint m_coefR;
  CFuint m_coefC;
  CFuint m_coefr;
  CFuint m_coefc;
  
  /// pointer to the values in this block
  CFreal* m_values;
  
}; // end of class BlockAccumulatorBaseCUDA
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BlockAccumulatorBaseCUDA_hh
