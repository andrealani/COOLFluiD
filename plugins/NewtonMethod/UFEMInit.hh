// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_U.minit_hh
#define COOLFluiD_Numerics_NewtonMethod_U.minit_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class U.minit : public NewtonIteratorCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);


  /// Constructor.
  explicit U.minit(const std::string& name);

  /// Destructor.
  ~U.minit()
  {
  }

  /// Execute Processing actions
  void execute();


  /// Configures this command with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private:

  /// Compute the value of the time limiter
  void computeTimeLimiter();


protected: // data

  /// storage of the past time rhs
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of the past time rhs
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// storage of the past time rhs
  Framework::DataSocketSink< Framework::State*> socket_pastPastStates;

  /// storage of the (time) limiter
  Framework::DataSocketSource< CFreal> socket_timeLimiter;

  /// time limiter
  Common::SelfRegistPtr<FiniteVolume::TimeLimiter> _timeLimiter;

  std::string _timeLimiterStr;

}; // class U.minit

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_U.minit_hh

