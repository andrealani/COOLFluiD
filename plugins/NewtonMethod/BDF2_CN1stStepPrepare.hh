// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_BDF2_CN1stStepPrepare_hh
#define COOLFluiD_Numerics_NewtonMethod_BDF2_CN1stStepPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This command computes the steady rhs of the initial solution,
  /// used with the BDF2 method, where for the first step Crank-Nicholson is used.
class BDF2_CN1stStepPrepare : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit BDF2_CN1stStepPrepare(std::string name) :
    NewtonIteratorCom(name),
    socket_rhs("rhs"),
    socket_pastTimeRhs("pastTimeRhs")
  {
  }

  /// Destructor.
  ~BDF2_CN1stStepPrepare()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// storage of the rhs
  Framework::DataSocketSink< CFreal> socket_rhs;

  /// storage of the past time rhs
  Framework::DataSocketSink< CFreal> socket_pastTimeRhs;

}; // class BDF2_CN1stStepPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_BDF2_CN1stStepPrepare_hh

