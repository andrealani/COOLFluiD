// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_StdParSolveSys_hh
#define COOLFluiD_Numerics_Trilinos_StdParSolveSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "Trilinos.hh"
#include "Framework/Storage.hh"
#include "TrilinosLSSData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a standard command to be executed during
   * solution of the linear system using Trilinos library
   */

//////////////////////////////////////////////////////////////////////////////

class StdParSolveSys : public TrilinosLSSCom {
public:

  /**
   * Constructor.
   */
  explicit StdParSolveSys(const std::string& name);

  /**
   * Destructor.
   */
  ~StdParSolveSys();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private: // data

  // socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  // socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;
  
  // local IDs of the locally owned unknowns (not the states !!)
  int *_lid;
  
  // global IDs of the locally owned unknowns (not the states !!)
  int *_gid;

}; // class SolveSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_StdParSolveSys_hh
