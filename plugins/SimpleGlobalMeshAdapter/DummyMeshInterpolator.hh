// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_DummyMeshInterpolator_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_DummyMeshInterpolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to read a newly created mesh
   *
   * @author Thomas Wuilbaut
   *
   */
class DummyMeshInterpolator : public SimpleMeshAdapterCom {
public:

  /**
   * Constructor.
   */
  explicit DummyMeshInterpolator(const std::string& name);

  /**
   * Destructor.
   */
  ~DummyMeshInterpolator()
  {
  }

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: // data

  /// Socket for states
  Framework::DataSocketSink <Framework::State* , Framework::GLOBAL> socket_states;

  /// Socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_otherStates;

}; // class DummyMeshInterpolator

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_DummyMeshInterpolator_hh

