// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CellVertexSparsityNoBlock_hh
#define COOLFluiD_Framework_CellVertexSparsityNoBlock_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalJacobianSparsity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Computes the sparticity of the global jacobian matrix for a cell-vertex
/// space discretization methods.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API CellVertexSparsityNoBlock : public GlobalJacobianSparsity {

public: // functions

  /// Constructor.
  CellVertexSparsityNoBlock();

  /// Destructor.
  virtual ~CellVertexSparsityNoBlock();

  /// Computes the non zero entries in the global jacobian matrix
  virtual void computeNNz(std::valarray<CFint>& nnz,
                          std::valarray<CFint>& ghostNnz);

  /// Computes the non zero entries in the global jacobian matrix.
  /// Also computes the vertex-vertex connectivity.
  virtual void computeMatrixPattern(
    std::valarray<CFint>& nnz,
    std::valarray<CFint>& ghostNnz,
    std::vector<std::vector<CFuint> >& matrixPattern);

private:

  /// computes the info on periodic pairs
  /// crossNodes: -1 is not periodic, otherwise the local idx of its pair
  void initPeriodicData(std::vector<int>& crossNodes);

}; // class CellVertexSparsityNoBlock

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CellVertexSparsityNoBlock_hh
