// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_NewmarkPrepare_hh
#define COOLFluiD_Numerics_NewtonMethod_NewmarkPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class NewmarkPrepare : public NewtonIteratorCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);


  /// Constructor.
  explicit NewmarkPrepare(std::string name);

  /// Destructor.
  ~NewmarkPrepare()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // handle to past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  // handle to past states
  Framework::DataSocketSink<Framework::State*> socket_pastStatesD;

  // handle to past states
  Framework::DataSocketSink<Framework::State*> socket_pastStatesD2;

  /// convergence parameter
  CFreal _alpha;

  /// convergence parameter
  CFreal _gamma;


}; // class NewmarkPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_NewmarkPrepare_hh

