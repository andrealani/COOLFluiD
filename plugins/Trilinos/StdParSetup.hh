// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_StdParSetup_hh
#define COOLFluiD_Numerics_Trilinos_StdParSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Trilinos.hh"
#include "TrilinosLSSData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

/**
 * Command to setup the Trilinos linear system solver Method
 * in a parallel run.
 */
class StdParSetup : public TrilinosLSSCom {
public:

  /**
   * Constructor.
   */
  explicit StdParSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdParSetup();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private: //helper functions

  /**
   *  Sets the Epetra communicator
   */
  void setEpetraComm();

  /**
   * Set up the index mapping
   */
  void setMapping();

  /**
   * Set up the matrix
   */
  void setMatrix();

  /**
   * Set up the vectors
   */
  void setVectors();

protected:

  /// socket for bStatesNeighbors
  /// It is a list of the neighbor states for the boundary states.
  /// It will be useful to avoid very expensive jacobian matrix
  /// reallocations when applying strong boundary condition
  Framework::DataSocketSink<std::valarray<Framework::State*> > socket_bStatesNeighbors;

  // socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  // socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_StdParSetup_hh

