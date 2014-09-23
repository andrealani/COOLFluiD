// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>
#include "Framework/BlockAccumulatorBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

BlockAccumulatorBase::BlockAccumulatorBase(const CFuint nbRows,
					   const CFuint nbCols,
					   const CFuint subBlockSize,
					   CFreal* array)
{ 
  reset(nbRows, nbCols, subBlockSize, array);
}

//////////////////////////////////////////////////////////////////////////////

BlockAccumulatorBase::~BlockAccumulatorBase()
{
}

//////////////////////////////////////////////////////////////////////////////
 
void BlockAccumulatorBase::reset(const CFuint nbRows,
				 const CFuint nbCols,
				 const CFuint subBlockSize,
				 CFreal* array)
{
  m_m = nbRows;
  m_n = nbCols;
  m_nb = subBlockSize;
  m_bsize = m_m*m_n*m_nb*m_nb;
  
  // assign the array pointer
  m_values = array;
  
  // default = row by row storage
  setCoefIndex(m_nb*m_nb*m_n, m_nb, m_nb*m_n, 1);
}
 
//////////////////////////////////////////////////////////////////////////////
    
void BlockAccumulatorBase::print() const
{
  // loop over row index of sub-blocks
  for (CFuint i = 0; i < m_m; ++i) {
    // loop over column index of sub-blocks
    for (CFuint j = 0; j < m_n; ++j) {
      std::cout << "(" << i << "," << j << ")"<< std::endl;
      for (CFuint ib = 0; ib < m_nb; ++ib) {
	for (CFuint jb = 0; jb < m_nb; ++jb) {
	  const CFuint idx = getIndex(i,j,ib,jb);
	  cf_assert(idx < m_bsize);
	  std::cout.precision(12);  std::cout << m_values[idx] << " ";
	}
	std::cout << std::endl;
      }
    }
  }
  std::cout << std::endl;
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
