// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GlobalJacobianSparsity_hh
#define COOLFluiD_Framework_GlobalJacobianSparsity_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include "Common/StringOps.hh"
#include "Common/OwnedObject.hh"
#include "Common/SafePtr.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {
    template <class T> class ConnectivityTable;
  }

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Interface of the classes that compute the sparticity of the
/// global jacobian matrix
/// Ususally used by LinearSystemSolver's to allocate correctly the global matrix,
/// possibly in parallel.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API GlobalJacobianSparsity : public Common::OwnedObject {

public: // functions

  typedef Environment::ConcreteProvider<GlobalJacobianSparsity> PROVIDER;

  /// Constructor.
  GlobalJacobianSparsity();

  /// Destructor.
  virtual ~GlobalJacobianSparsity();

  /// Computes the non zero entries in the global jacobian matrix
  virtual void computeNNz(std::valarray<CFint>& nnz,
                          std::valarray<CFint>& ghostNnz) = 0;
  
  /// Computes the non zero entries in the global jacobian matrix
  /// using nodes instead of states
  virtual void computeNNzNodeBased(std::valarray<CFint>& nnz,
				   std::valarray<CFint>& ghostNnz) 
  {
    throw Common::NotImplementedException
      (FromHere(), "GlobalJacobianSparsity::computeNNzNodeBased()");
  }
  
  /// Computes the non zero entries in the global jacobian matrix.
  /// Also computes the vertex-vertex connectivity.
  virtual void computeMatrixPattern(
    std::valarray<CFint>& nnz,
    std::valarray<CFint>& ghostNnz,
    std::vector<std::vector<CFuint> >& matrixPattern) = 0;

  /// Computes the matrix patern and stores it in a connectivity table
  virtual void computeMatrixPattern
  (DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
   Common::ConnectivityTable<CFuint>& matrixPattern);

  /// Gets the Class name
  static std::string getClassName() { return "GlobalJacobianSparsity";  }

  /// Sets the DataSockets
  void setDataSockets
  (DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
   DataSocketSink<Framework::Node*, Framework::GLOBAL> nodesSocket,
   DataSocketSink<std::valarray<State*> > bStatesNeighborsSocket);

protected: // helper function

  /// Set the array of flags telling if a State is on the boundary
  void computeBoundaryStatesFlag(std::valarray<bool>& isBoundaryState) const;

protected: // data
  
  /// socket for States
  DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for Nodes
  DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// socket for boundaryStates neighbors
  DataSocketSink<std::valarray<State*> > socket_bStatesNeighbors;
  
}; // class GlobalJacobianSparsity

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(GlobalJacobianSparsity)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GlobalJacobianSparsity_hh

