// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Pardiso_StdSetup_hh
#define COOLFluiD_Pardiso_StdSetup_hh

#include "Pardiso/PardisoData.hh"
#include "Framework/DataSocketSink.hh"

namespace COOLFluiD {
  namespace Pardiso {

/// This is a standard command to setup the Pardiso method
class StdSetup : public PardisoCom {

public:

  /// Constructor
  explicit StdSetup(const std::string& name) :
    PardisoCom(name),
    socket_bStatesNeighbors("bStatesNeighbors"),
    socket_states("states") {}

  /// Destructor
  ~StdSetup() {}

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSockets that this command needs as sinks
   * @return vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets() {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > result;
    result.push_back(&socket_states);
    result.push_back(&socket_bStatesNeighbors);
    return result;
  }


protected:

  /**
   * socket for bStatesNeighbors
   * It's a list of neighbor states for the boundary states (to avoid matrix
   * reallocations when applying boundary conditions)
   */
  Framework::DataSocketSink<
    std::valarray<Framework::State*> > socket_bStatesNeighbors;

  // socket for states
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL >
    socket_states;


private: // methods

  /// Get state-state connectivity to build matrix structure
  void getStructure(
    std::vector< CFint >& nnz,
    std::vector< std::vector< CFuint > >& nz );

  /// Expand non-zeros and number of non-zeros mapping
  void expandStructureMapping(
    std::vector< std::vector< CFuint > >& nz1, std::vector< CFint >& nnz1,
    const CFuint nbEqs );

}; // class StdSetup


  } // namespace Pardiso
} // namespace COOLFluiD

#endif // COOLFluiD_Numerics_Pardiso_StdSetup_hh

