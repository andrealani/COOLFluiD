// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_StdParSolveSys_hh
#define COOLFluiD_Numerics_Petsc_StdParSolveSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"

#include "Petsc/PetscLSSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class State; }

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a standard command to be executed during
  * solution of the linear system using Petsc library
  *
  * @author Andrea Lani
  *
  */
class StdParSolveSys : public PetscLSSCom {
public:

  /**
   * Constructor.
   */
  explicit StdParSolveSys(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdParSolveSys();

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected: // data
  
  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// socket for the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;
  
  /// vector of the local ids of only the updatable states
  std::vector<CFint> _upLocalIDs;
  
  /// indexes for the insertion of elements in a PetscVector
  std::vector<CFint> _upStatesGlobalIDs;
   
}; // class SolveSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_StdParSolveSys_hh
