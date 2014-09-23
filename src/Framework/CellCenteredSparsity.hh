// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CellCenteredSparsity_hh
#define COOLFluiD_Framework_CellCenteredSparsity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalJacobianSparsity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Computes the sparticity of the global jacobian matrix for a cell-centered
/// space discretization methods.
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Kris Van den Abeele
class Framework_API CellCenteredSparsity : public GlobalJacobianSparsity {

public: // functions

  /// Constructor.
  CellCenteredSparsity();

  /// Destructor.
  virtual ~CellCenteredSparsity();

  /// Computes the non zero entries in the global jacobian matrix
  virtual void computeNNz(std::valarray<CFint>& nnz,
                          std::valarray<CFint>& ghostNnz);

  /// Computes the non zero entries in the global jacobian matrix.
  /// Also computes the vertex-vertex connectivity.
  virtual void computeMatrixPattern(
    std::valarray<CFint>& nnz,
    std::valarray<CFint>& ghostNnz,
    std::vector<std::vector<CFuint> >& matrixPattern);
  
  /// Computes the matrix patern and stores it in a connectivity table
  virtual void computeMatrixPattern
  (DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
   Common::ConnectivityTable<CFuint>& matrixPattern);
  
}; // class CellCenteredSparsity

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CellCenteredSparsity_hh
