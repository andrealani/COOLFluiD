// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_BaseSetup_hh
#define COOLFluiD_Numerics_Petsc_BaseSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header
#include "Petsc/PetscLSSData.hh"

#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class is the a base class for commands that setup Petsc
  * Method. The execute() function is a template method. Subclasses
  * are meant to simply override the (private) hook methods
  *
  * @author Tiago Quintino
  * @author Andrea Lani
  */
class BaseSetup : public PetscLSSCom {
public:

  /**
   * Constructor.
   */
  explicit BaseSetup(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~BaseSetup();

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /**
   * Set up the KSP solver
   */
  virtual void setKSP();

  /**
   * Set up the index mapping
   */
  virtual void setIdxMapping() = 0;

  /**
   * Set up the matrix
   */
  virtual void setMatrix(const CFuint localSize,
			 const CFuint globalSize) = 0;

  /**
   * Set up the vectors
   */
  virtual void setVectors(const CFuint localSize,
			  const CFuint globalSize);

protected: // data

  /// socket for bStatesNeighbors
  /// It is a list of the neighbor states for the boundary states.
  /// It will be useful to avoid very expensive jacobian matrix
  /// reallocations when applying strong boundary condition
  Framework::DataSocketSink<std::valarray<Framework::State*> > socket_bStatesNeighbors;

  /// socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// socket for the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
   /// socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_BaseSetup_hh

