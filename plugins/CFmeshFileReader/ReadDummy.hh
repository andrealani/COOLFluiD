// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_ReadDummy_hh
#define COOLFluiD_CFmeshFileReader_ReadDummy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"

#include "CFmeshFileReader/CFmeshReaderData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class is a MethodCommand implementing a dummy reader
/// @author Andrea Lani
class CFmeshFileReader_API ReadDummy : public CFmeshReaderCom {

public: // functions
  
  /// Constructor.
  explicit ReadDummy(const std::string& name);

  /// Destructor.
  virtual ~ReadDummy();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /// Configures the command.
  void configure ( Config::ConfigArgs& args );
  
  /// Execute Processing actions
  void execute();
  
protected: // data
 
  /// The set of data sockets to be used by the strategy
  Framework::DynamicDataSocketSet<> m_sockets;
  
  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*,Framework::GLOBAL> socket_nodes;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*,Framework::GLOBAL> socket_states;

}; // class ReadDummy

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileReader_ReadDummy_hh
