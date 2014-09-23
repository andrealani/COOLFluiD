// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_CopySol_hh
#define COOLFluiD_Numerics_ForwardEuler_CopySol_hh

//////////////////////////////////////////////////////////////////////////////

#include "FwdEulerData.hh"
#include "Framework/DataSocketSink.hh"

namespace COOLFluiD {

  namespace Framework {
    class State;
  }


    namespace ForwardEuler {

/// This command copies the right-hand side vector to the States
class  ForwardEuler_API CopySol : public FwdEulerCom {
public:

  /// Constructor
  explicit CopySol(const std::string& name) :
    FwdEulerCom(name),
    socket_states("states"),
    socket_rhs("rhs")
  {
  }

  /// Destructor
  ~CopySol()
  {
  }

  /// Execute processing actions
  void execute();

  /// @return a vector of SafePtr with the DataSockets needed as sinks
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();


private:

  /// Handle to States
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// Handle to RHS
  Framework::DataSocketSink<CFreal >   socket_rhs;

};



  }  // namespace Numerics
}  // namespace COOLFluiD

#endif // COOLFluiD_Numerics_ForwardEuler_CopySol_hh

