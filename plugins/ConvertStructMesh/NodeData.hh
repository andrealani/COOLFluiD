// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshTools_NodeData_hh
#define COOLFluiD_CFmeshTools_NodeData_hh

//////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <fstream>

#include "MathTools/RealMatrix.hh"
#include "Common/CFMap3D.hh"
#include "Common/SafePtr.hh"
#include "ConvertStructMesh/Block.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers and interface to handle node data
 *
 * @author Andrea Lani
 *
 */
typedef Common::CFMap3D<CFreal, CFreal, CFreal, CFuint*> MapCoord2ID;

class Node3D {
public:
  
  Node3D(CFreal* x) : _x(x) {}

  void reset(CFreal* x)
  {
    cf_assert(x != CFNULL);
    _x = x;
  }

  friend bool operator< (const Node3D& n1, const Node3D& n2)
  {
    cf_assert(n1._x != CFNULL);
    cf_assert(n2._x != CFNULL);
    if (n1._x[XX] < n2._x[XX]) return true;
    if (n1._x[XX] > n2._x[XX]) return false;
    if (n1._x[YY] < n2._x[YY]) return true;
    if (n1._x[YY] > n2._x[YY]) return false;
    return (n1._x[ZZ] < n2._x[ZZ]);
  }

  friend bool operator> (const Node3D& n1, const Node3D& n2)
  {
    cf_assert(n1._x != CFNULL);
    cf_assert(n2._x != CFNULL);
    if (n1._x[XX] > n2._x[XX]) return true;
    if (n1._x[XX] < n2._x[XX]) return false;
    if (n1._x[YY] > n2._x[YY]) return true;
    if (n1._x[YY] < n2._x[YY]) return false;
    return (n1._x[ZZ] > n2._x[ZZ]);
  }

  friend bool operator!= (const Node3D& n1, const Node3D& n2)
  {
    return ((n1._x[XX] != n2._x[XX]) ||
	    (n1._x[YY] != n2._x[YY]) ||
	    (n1._x[ZZ] != n2._x[ZZ]));
  }

  friend bool operator== (const Node3D& n1, const Node3D& n2)
  {
    return !(operator!=(n1,n2));
  }
  
  /**
   * Overloading of the stream operator "<<" for the output
   */
  friend std::ostream& operator<<(std::ostream& out, const Node3D& n)
  {
    cf_assert(n._x != CFNULL);
    using namespace std;

    out << setw(16) << fixed << setprecision(8)
	<< n._x[XX] << " " << n._x[YY] << " " << n._x[ZZ] << "\n";
    return out;
  }

  /**
   * Overloading of the stream operator[]
   */
  CFreal operator[] (CFuint iDim) const
  {
    cf_assert(_x != CFNULL);
    return _x[iDim];
  }

private:
  CFreal* _x;
};

class NodeData {
public:
  /**
   * Constructor
   */
  NodeData() : _dim(0), _nx(), _ny(), _nz(), _allNodes(), _nodeIDs()
  {
  }

  /**
   * Default destructor
   */
  ~NodeData()
  {
  }

  /**
   * Set the space dimension
   */
  void setDim(CFuint dim)
  {
    cf_assert(dim > 1);
    _dim = dim;
  }

  /**
   * Set the number of blocks
   */
  void setNbBlocks(CFuint nbBlocks)
  {
    cf_assert(nbBlocks > 0);
    _nodeIDs.resize(nbBlocks);
  }

  /**
   * Get the global ID of the node having indices (i,j,k) in block ib
   */
  CFuint getNodeID(CFuint ib, CFuint i, CFuint j, CFuint k) const
  {
    cf_assert(ib < _nodeIDs.size());
    return _nodeIDs[ib][k*_nx[ib]*_ny[ib] + j*_nx[ib] + i];
  }

  /**
   * Get the global ID of the node having indices (i,j,k) in block ib
   */
  CFuint getNodeID(CFuint ib, CFuint i, CFuint j) const
  {
    cf_assert(ib < _nodeIDs.size());
    return _nodeIDs[ib][j*_nx[ib] + i];
  }

  /**
   * Get the array storing the maximum number of block nodes in x direction
   */
  std::vector<CFuint>& getNxVec() {return _nx;}

  /**
   * Get the array storing the maximum number of block nodes in y direction
   */
  std::vector<CFuint>& getNyVec() {return _ny;}

  /**
   * Get the array storing the maximum number of block nodes in z direction
   */
  std::vector<CFuint>& getNzVec() {return _nz;}

  /**
   * Assign the node IDs
   */
  void assignIDs(MapCoord2ID& bcoord2LocalID,
	         const std::vector<std::vector<bool> >& isBnode,
		 std::vector<RealMatrix>& xyzInBlock,
		 CFuint nbBNodes);
  /**
   * Assign the node IDs
   */
  void assignIDs(std::vector<Block>& blocks,
	         const std::vector<std::vector<bool> >& isBnode,
		 std::vector<RealMatrix>& xyzInBlock);

  /**
   * Assign defualt IDs
   */
  void assignDefaultIDs(CFuint iBlock, CFuint nbNodes);

  /**
   * Get the full node storage
   */
  Common::SafePtr<RealMatrix> getData()
  {
    return &_allNodes;
  }

private: // data

  /// space dimension
  CFuint _dim;

  /// array storing the maximum number of block nodes in x direction
  std::vector<CFuint> _nx;

  /// array storing the maximum number of block nodes in y direction
  std::vector<CFuint> _ny;

  /// array storing the maximum number of block nodes in z direction
  std::vector<CFuint> _nz;

  /// all nodes
  RealMatrix _allNodes;

  /// node IDs per block
  std::vector<std::vector<int> > _nodeIDs;

}; // end of class NodeData

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshTools_NodeData_hh
