// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_Tecplot2CFmesh_NodeDim_hh
#define COOLFluiD_IO_Tecplot2CFmesh_NodeDim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "MathTools/MathChecks.hh"
#include <cassert>
#include <iomanip>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Tecplot2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers and interface to handle node data
 *
 * @author Andrea Lani
 *
 */
 
//////////////////////////////////////////////////////////////////////////////
     
class NodeDim {
public:
  
  /// Constructor
  NodeDim(CFreal* x) : m_x(x) {}
  
  /// Reset the pointer
  void reset(CFreal* x)
  {
    assert(x != CFNULL);
    m_x = x;
  }
  
  /// Overloading of the operator < for comparison
  friend bool operator< (const NodeDim& n1, const NodeDim& n2)
  {
    assert(n1.m_x != CFNULL);
    assert(n2.m_x != CFNULL);
    if (n1.m_x[XX] < n2.m_x[XX]) return true;
    if (n1.m_x[XX] > n2.m_x[XX]) return false;
    if (m_ssize == 2) return (n1.m_x[YY] < n2.m_x[YY]);
    if (n1.m_x[YY] < n2.m_x[YY]) return true;
    if (n1.m_x[YY] > n2.m_x[YY]) return false;
    return (n1.m_x[ZZ] < n2.m_x[ZZ]);
  }
  
  /// Overloading of the operator > for comparison
  friend bool operator> (const NodeDim& n1, const NodeDim& n2)
  {
    assert(n1.m_x != CFNULL);
    assert(n2.m_x != CFNULL);
    if (n1.m_x[XX] > n2.m_x[XX]) return true;
    if (n1.m_x[XX] < n2.m_x[XX]) return false;
    if (m_ssize == 2) return (n1.m_x[YY] > n2.m_x[YY]);
    if (n1.m_x[YY] > n2.m_x[YY]) return true;
    if (n1.m_x[YY] < n2.m_x[YY]) return false;
    return (n1.m_x[ZZ] > n2.m_x[ZZ]);
  }
  
  /// Overloading of the operator != for comparison
  friend bool operator!= (const NodeDim& n1, const NodeDim& n2)
  {
    if (m_ssize == 3) {
      // return ((n1.m_x[XX] != n2.m_x[XX]) ||
      // 	       (n1.m_x[YY] != n2.m_x[YY]) ||
      // 	       (n1.m_x[ZZ] != n2.m_x[ZZ]));
      return (MathTools::MathChecks::isNotEqual(n1.m_x[XX], n2.m_x[XX]) ||
	      MathTools::MathChecks::isNotEqual(n1.m_x[YY], n2.m_x[YY]) ||
	      MathTools::MathChecks::isNotEqual(n1.m_x[ZZ], n2.m_x[ZZ]));
    }
    // return ((n1.m_x[XX] != n2.m_x[XX]) || (n1.m_x[YY] != n2.m_x[YY]));
    return (MathTools::MathChecks::isNotEqual(n1.m_x[XX], n2.m_x[XX]) || 
	    MathTools::MathChecks::isNotEqual(n1.m_x[YY], n2.m_x[YY]));
  }
  
  /// Overloading of the operator == for comparison
  friend bool operator== (const NodeDim& n1, const NodeDim& n2)
  {
    return !(operator!=(n1,n2));
  }
  
  /// Overloading of the stream operator "<<" for the output
  friend std::ostream& operator<<(std::ostream& out, const NodeDim& n)
  {
    assert(n.m_x != CFNULL);
    using namespace std;
    for (CFuint iDim = 0; iDim < m_ssize; ++iDim) {
      out << setw(22) << fixed << setprecision(14) << n.m_x[iDim] << " ";
    }
    out << "\n";
    return out;
  }
  
  
  /// Overloading of the stream operator[]
  CFreal operator[] (CFuint iDim) const
  {
    assert(m_x != CFNULL);
    return m_x[iDim];
  }
  
  /// Overloading of the stream operator[]
  CFreal& operator[] (CFuint iDim)
  {
    assert(m_x != CFNULL);
    return m_x[iDim];
  }
  
  /// @return the size of this array
  static CFuint size() 
  {
    cf_assert(m_ssize > 0); 
    return m_ssize;
  }
  
  /// set the size of this array
  static void setSize(const CFuint ssize) 
  {
    cf_assert(ssize > 0);  
    m_ssize = ssize;
  }
  
private:
  /// data pointer
  CFreal* m_x;
  
  /// static size to be shared by all @see NodeDim
  static CFuint  m_ssize;
  
}; // end class NodeDim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Tecplot2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Tecplot2CFmesh_NodeDim_hh
