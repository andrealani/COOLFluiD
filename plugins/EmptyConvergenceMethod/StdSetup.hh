// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptyConvergenceMethod_StdSetup_hh
#define COOLFluiD_EmptyConvergenceMethod_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

#include "EmptyConvergenceMethod/EmptyIteratorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class State;
  }
  
  namespace EmptyConvergenceMethod {
    
//////////////////////////////////////////////////////////////////////////////

/// This class represents a standard setup for an  empty ConvergenceMethod
/// @author Andrea Lani
class  EmptyConvergenceMethod_API StdSetup : public EmptyIteratorCom {
public:

  /// Constructor.
  explicit StdSetup(const std::string& name);

  /// Destructor.
  virtual ~StdSetup();
  
  /// Configures this Command
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /// Execute Processing actions
  virtual void execute();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for Rhs
  Framework::DataSocketSource<CFreal> socket_rhs;

  /// socket for updateCoeff, denominators of the coefficients for the update
  Framework::DataSocketSource<CFreal> socket_updateCoeff;
  
  /// socket of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace EmptyConvergenceMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_EmptyConvergenceMethod_StdSetup_hh

