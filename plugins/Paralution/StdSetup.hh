// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Paralution_StdSetup_hh
#define COOLFluiD_Numerics_Paralution_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Paralution/ParalutionLSSData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class is the a base class for commands that setup Paralution
  * Method. The execute() function is a template method. Subclasses
  * are meant to simply override the (private) hook methods
  *
  * @author Isaac Alonso
  * @author Andrea Lani
  */
class StdSetup : public ParalutionLSSCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdSetup();

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
   * Set up the Krylov solver
   */
  virtual void setKSP();
  
  /**
   * Set up the index mapping
   */
  virtual void setIdxMapping();
  
  /**
   * Set up the matrix
   */
  virtual void setMatrix(const CFuint localSize,
			 const CFuint globalSize);

  /**
   * Set up the vectors
   */
  virtual void setVectors(const CFuint localSize,
			  const CFuint globalSize);
protected: // data
  
  /// socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// socket for the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Paralution_StdSetup_hh

