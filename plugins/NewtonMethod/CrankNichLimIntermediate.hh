// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_CrankNichLimIntermediate_hh
#define COOLFluiD_Numerics_NewtonMethod_CrankNichLimIntermediate_hh

//////////////////////////////////////////////////////////////////////////////

#include "CrankNichIntermediate.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
  /// @author Thomas Wuilbaut
class CrankNichLimIntermediate : public CrankNichIntermediate {
public:

  /// Constructor.
  explicit CrankNichLimIntermediate(std::string name) :
    CrankNichIntermediate(name),
    socket_timeLimiter("timeLimiter")

  {
  }

  /// Destructor.
  ~CrankNichLimIntermediate()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  // socket for the time limiter
  Framework::DataSocketSink<CFreal> socket_timeLimiter;

}; // class CrankNichLimIntermediate

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_CrankNichLimIntermediate_hh

