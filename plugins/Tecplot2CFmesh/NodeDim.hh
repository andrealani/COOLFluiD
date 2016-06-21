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
template <CFuint DIM> class NodeDim{};
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This explicit specialization offers and interface to handle 3D node data
 *
 * @author Andrea Lani
 *
 */
template <> 
class NodeDim<3> {
public:
  
  /// Constructor
  NodeDim(CFreal* x) : _x(x) {}
  
  /// Reset the pointer
  void reset(CFreal* x)
  {
    assert(x != CFNULL);
    _x = x;
  }
  
  /// Overloading of the operator < for comparison
  friend bool operator< (const NodeDim<3>& n1, const NodeDim<3>& n2)
  {
    assert(n1._x != CFNULL);
    assert(n2._x != CFNULL);
    if (n1._x[XX] < n2._x[XX]) return true;
    if (n1._x[XX] > n2._x[XX]) return false;
    if (n1._x[YY] < n2._x[YY]) return true;
    if (n1._x[YY] > n2._x[YY]) return false;
    return (n1._x[ZZ] < n2._x[ZZ]);
  }
  
  /// Overloading of the operator > for comparison
  friend bool operator> (const NodeDim<3>& n1, const NodeDim<3>& n2)
  {
    assert(n1._x != CFNULL);
    assert(n2._x != CFNULL);
    if (n1._x[XX] > n2._x[XX]) return true;
    if (n1._x[XX] < n2._x[XX]) return false;
    if (n1._x[YY] > n2._x[YY]) return true;
    if (n1._x[YY] < n2._x[YY]) return false;
    return (n1._x[ZZ] > n2._x[ZZ]);
  }
  
  /// Overloading of the operator != for comparison
  friend bool operator!= (const NodeDim<3>& n1, const NodeDim<3>& n2)
  {
    // return ((n1._x[XX] != n2._x[XX]) ||
    // 	    (n1._x[YY] != n2._x[YY]) ||
    // 	    (n1._x[ZZ] != n2._x[ZZ]));
    
    return (MathTools::MathChecks::isNotEqual(n1._x[XX], n2._x[XX]) ||
	    MathTools::MathChecks::isNotEqual(n1._x[YY], n2._x[YY]) ||
	    MathTools::MathChecks::isNotEqual(n1._x[ZZ], n2._x[ZZ]));
  }
  
  /// Overloading of the operator == for comparison
  friend bool operator== (const NodeDim<3>& n1, const NodeDim<3>& n2)
  {
    return !(operator!=(n1,n2));
  }
  
  /// Overloading of the stream operator "<<" for the output
  friend std::ostream& operator<<(std::ostream& out, const NodeDim<3>& n)
  {
    assert(n._x != CFNULL);
    using namespace std;
    out << setw(22) << fixed << setprecision(14)
	<< n._x[XX] << " " << n._x[YY] << " " << n._x[ZZ] << "\n";
    return out;
  }
  
  
  /// Overloading of the stream operator[]
  CFreal operator[] (CFuint iDim) const
  {
    assert(_x != CFNULL);
    return _x[iDim];
  }
  
  /// Overloading of the stream operator[]
  CFreal& operator[] (CFuint iDim)
  {
    assert(_x != CFNULL);
    return _x[iDim];
  }
  
  /// @return the size of this array
  CFuint size() const {return 3;}
  
private:
  CFreal* _x;
  
}; // end class NodeDim

//////////////////////////////////////////////////////////////////////////////

/**
 * This explicit specialization offers and interface to handle 2D node data
 *
 * @author Andrea Lani
 *
 */
template <> 
class NodeDim<2> {
public:
  
  /// Constructor
  NodeDim(CFreal* x) : _x(x) {}
  
  /// Reset the pointer
  void reset(CFreal* x)
  {
    assert(x != CFNULL);
    _x = x;
  }
  
  /// Overloading of the operator < for comparison
  friend bool operator< (const NodeDim<2>& n1, const NodeDim<2>& n2)
  {
    assert(n1._x != CFNULL);
    assert(n2._x != CFNULL);
    if (n1._x[XX] < n2._x[XX]) return true;
    if (n1._x[XX] > n2._x[XX]) return false;
    return (n1._x[YY] < n2._x[YY]);
  }
  
  /// Overloading of the operator > for comparison
  friend bool operator> (const NodeDim<2>& n1, const NodeDim<2>& n2)
  {
    assert(n1._x != CFNULL);
    assert(n2._x != CFNULL);
    if (n1._x[XX] > n2._x[XX]) return true;
    if (n1._x[XX] < n2._x[XX]) return false;
    return (n1._x[YY] > n2._x[YY]);
  }
  
  /// Overloading of the operator != for comparison
  friend bool operator!= (const NodeDim<2>& n1, const NodeDim<2>& n2)
  {
    // return ((n1._x[XX] != n2._x[XX]) || (n1._x[YY] != n2._x[YY]));
    return (MathTools::MathChecks::isNotEqual(n1._x[XX], n2._x[XX]) || 
	    MathTools::MathChecks::isNotEqual(n1._x[YY], n2._x[YY]));
  }
  
  /// Overloading of the operator == for comparison
  friend bool operator== (const NodeDim<2>& n1, const NodeDim<2>& n2)
  {
    return !(operator!=(n1,n2));
  }
  
  /// Overloading of the stream operator "<<" for the output
  friend std::ostream& operator<<(std::ostream& out, const NodeDim<2>& n)
  {
    assert(n._x != CFNULL);
    using namespace std;
    out << setw(22) << fixed << setprecision(14) 
	<< n._x[XX] << " " << n._x[YY] << "\n";
    return out;
  }
  
  
  /// Overloading of the stream operator[]
  CFreal operator[] (CFuint iDim) const 
  {    
    assert(_x != CFNULL);
    return _x[iDim];
  }
  
  /// Overloading of the stream operator[]
  CFreal& operator[] (CFuint iDim)
  {
    assert(_x != CFNULL);
    return _x[iDim];
  }
  
  /// @return the size of this array
  CFuint size() const {return 2;}
  
private:
  CFreal* _x;

}; // end class NodeDim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Tecplot2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Tecplot2CFmesh_NodeDim_hh
