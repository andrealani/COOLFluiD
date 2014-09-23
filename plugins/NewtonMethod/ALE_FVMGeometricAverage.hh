// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_ALE_FVMGeometricAverage_hh
#define COOLFluiD_Numerics_NewtonMethod_ALE_FVMGeometricAverage_hh

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
class ALE_FVMGeometricAverage : public NewtonIteratorCom {
public:

  /// Constructor.
  explicit ALE_FVMGeometricAverage(const std::string& name) :
    NewtonIteratorCom(name),
    socket_nodes("nodes"),
    socket_futureNodes("futureNodes")
  {
  }

  /// Destructor.
  ~ALE_FVMGeometricAverage()
  {
  }

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  // handle to nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // handle to futureNodes
  Framework::DataSocketSink<Framework::Node*> socket_futureNodes;


}; // class ALE_FVMGeometricAverage

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_ALE_FVMGeometricAverage_hh

