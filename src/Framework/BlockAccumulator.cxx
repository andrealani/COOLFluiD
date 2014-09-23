// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Framework/BlockAccumulator.hh"

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

BlockAccumulator::BlockAccumulator(const CFuint nbRows,
                                   const CFuint nbCols,
                                   const CFuint subBlockSize,
                                   const LSSIdxMapping& localToGlobal,
				   CFreal* ptr) :
  _localToGlobal(localToGlobal),
  m(nbRows),
  n(nbCols),
  nb(subBlockSize)
{
  m_bsize = m*n*nb*nb;
  
  if (m > 0) im.resize(m);
  if (n > 0) in.resize(n);
  
  if (ptr == CFNULL) {
    cf_assert(m_bsize > 0);
    values.resize(m_bsize, 0.);
  }
  else {
    values.wrap(m_bsize, ptr);
  }
  
  setCoefIndex(nb*nb*n, nb, nb*n, 1); // default = row by row storage
}
    
//////////////////////////////////////////////////////////////////////////////

BlockAccumulator::~BlockAccumulator()
{
}

//////////////////////////////////////////////////////////////////////////////

void BlockAccumulator::resetB(const CFuint nbRows,
                              const CFuint nbCols,
                              const CFuint subBlockSize,
                              CFreal* array)
{
  if (nbRows != m) {im.resize(nbRows);}
  if (nbCols != n) {in.resize(nbCols);}
  m = nbRows; 
  n = nbCols;
  nb = subBlockSize;
  values.wrap(m_bsize, array);
  setCoefIndex(nb*nb*n, nb, nb*n, 1); // default = row by row storage
}
    
//////////////////////////////////////////////////////////////////////////////

void BlockAccumulator::printToScreen() const
{
  CFout << "Block row indices   : ";
  for (CFuint i=0; i<m; i++) CFout << im[i] << " "; CFout << "\n";
  CFout << "Block col indices   : ";
  for (CFuint i=0; i<n; i++) CFout << in[i] << " "; CFout << "\n";
  for (CFuint R=0; R<m; R++) {
    for (CFuint r=0; r<nb; r++) {
      for (CFuint C=0; C<n; C++) {
        for (CFuint c=0; c<nb; c++) {
          CFout << " " << values[getIndex(R,C,r,c)];
        }
        CFout << " | ";
      }
      CFout << "\n";
    }
    CFout << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
