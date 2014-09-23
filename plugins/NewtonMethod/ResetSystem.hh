// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_ResetSystem_hh
#define COOLFluiD_Numerics_NewtonMethod_ResetSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand of the NewtonIterator method
  /// used for resetting the linear system matrix, and right hand side vector
  /// to zero.
class ResetSystem : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit ResetSystem(std::string name) :
    NewtonIteratorCom(name),
    socket_rhs("rhs")
  {
  }

  /// Destructor.
  ~ResetSystem()
  {
  }

  /// Execute Processing actions
  virtual void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  // socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

}; // class ResetSystem

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_ResetSystem_hh

