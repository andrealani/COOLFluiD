// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_BDF2_CN1stStepIntermediate_hh
#define COOLFluiD_Numerics_NewtonMethod_BDF2_CN1stStepIntermediate_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.

class BDF2_CN1stStepIntermediate : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit BDF2_CN1stStepIntermediate(std::string name) :
    NewtonIteratorCom(name),
    socket_rhs("rhs"),
    socket_pastTimeRhs("pastTimeRhs")
  {
  }

  /// Destructor.
  ~BDF2_CN1stStepIntermediate()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // variables

  /// storage of the rhs
  Framework::DataSocketSink< CFreal> socket_rhs;

  /// storage of the past time rhs
  Framework::DataSocketSink< CFreal> socket_pastTimeRhs;

}; // class BDF2_CN1stStepIntermediate

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_BDF2_CN1stStepIntermediate_hh

