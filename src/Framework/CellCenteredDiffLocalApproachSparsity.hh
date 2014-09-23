// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CellCenteredDiffLocalApproachSparsity_hh
#define COOLFluiD_Framework_CellCenteredDiffLocalApproachSparsity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalJacobianSparsity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Computes the sparticity of the global jacobian matrix for cell-centered
/// space discretization methods using a local approach for the diffusive terms.
/// @author Kris Van den Abeele
class Framework_API CellCenteredDiffLocalApproachSparsity : public GlobalJacobianSparsity {

public: // functions

  /// Constructor.
  CellCenteredDiffLocalApproachSparsity();

  /// Destructor.
  virtual ~CellCenteredDiffLocalApproachSparsity();

  /// Computes the non zero entries in the global jacobian matrix
  virtual void computeNNz(std::valarray<CFint>& nnz,
                          std::valarray<CFint>& ghostNnz);

  /// Computes the non zero entries in the global jacobian matrix.
  /// Also computes the vertex-vertex connectivity.
  virtual void computeMatrixPattern(
    std::valarray<CFint>& nnz,
    std::valarray<CFint>& ghostNnz,
    std::vector<std::vector<CFuint> >& matrixPattern);

}; // class CellCenteredDiffLocalApproachSparsity

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CellCenteredDiffLocalApproachSparsity_hh
