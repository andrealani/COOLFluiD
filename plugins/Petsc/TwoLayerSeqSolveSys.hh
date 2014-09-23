// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_TwoLayerSeqSolveSys_hh
#define COOLFluiD_Numerics_Petsc_TwoLayerSeqSolveSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Petsc/PetscLSSData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a standard command to be executed during
  * solution of the linear system using Petsc library
  *
  * @author Thomas Wuilbaut
  *
  */
class TwoLayerSeqSolveSys : public PetscLSSCom {
public:

  /**
   * Constructor.
   */
  explicit TwoLayerSeqSolveSys(const std::string& name) :
    PetscLSSCom(name),
    socket_states("states"),
    socket_rhs("rhs"),
    socket_interRhs("interRhs"),
    _idx()
  {
  }

  /**
   * Destructor.
   */
  ~TwoLayerSeqSolveSys()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * SetUp
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // handle to rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // handle to update coefficient
  Framework::DataSocketSink<CFreal> socket_interRhs;

  // array to store the local (=global) IDs
  std::vector<CFint> _idx;

}; // class SolveSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_TwoLayerSeqSolveSys_hh
