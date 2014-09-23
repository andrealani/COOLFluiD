// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_ReadBase_hh
#define COOLFluiD_CFmeshFileReader_ReadBase_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CFmeshReaderSource.hh"
#include "Framework/DataSocketSink.hh"

#include "CFmeshFileReader/CFmeshReaderData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class is a MethodCommand that reads the data in the file.
/// It serves as a base class for more specialized commands.
/// @author Tiago Quintino
class CFmeshFileReader_API ReadBase : public CFmeshReaderCom {

public: // functions


  /// Constructor.
  explicit ReadBase(const std::string& name);

  /// Destructor.
  virtual ~ReadBase();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets();

  /// Configures the command.
  void configure ( Config::ConfigArgs& args );

protected: // functions

  /// Corrects the size of the states in CFmeshData to the current
  /// PhysicalModelStack::getActive()->getNbEq()
  void correctStates();

  /// Applies the scalings to the nodes.
  /// One scaling will be permanent and is configurable by the user, affecting the output of the mesh.
  /// The other is based on the reference lenght for the simulation and is internal, thus it does
  /// not affect how the output mesh looks like. Output formatters have the responsability of scaling back.
  void applyScalings();

  /// Applies the a translation vector which is configurable to all the nodes.
  void applyTranslation();

protected: // data

  /// The set of data sockets to be used by the strategy
  Framework::DynamicDataSocketSet<> m_sockets;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*,Framework::GLOBAL> socket_nodes;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*,Framework::GLOBAL> socket_states;

  /// the actual data to write to the file
  std::auto_ptr<Framework::CFmeshReaderSource> m_data;

}; // class ReadBase

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileReader_ReadBase_hh
