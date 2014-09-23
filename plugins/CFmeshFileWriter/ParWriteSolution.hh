// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileWriter_ParWriteSolution_hh
#define COOLFluiD_CFmeshFileWriter_ParWriteSolution_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CFmeshWriterSource.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "CFmeshFileWriter/CFmeshWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// sent to Domain to be executed in order to setup the MeshData.
/// @author Andrea Lani
template <typename WRITER>
class CFmeshFileWriter_API ParWriteSolution : public CFmeshWriterCom {
public:

  /// Constructor.
  explicit ParWriteSolution(const std::string& name);

  /// Destructor.
  ~ParWriteSolution();

  /// Set up some private data needed
  /// for the simulation
  void setup();

  /// Configures the command.
  void configure ( Config::ConfigArgs& args );

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data

  // the actual data to write to the file
  std::auto_ptr<Framework::CFmeshWriterSource> _data;

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;
  
  /// CFmesh file writer
  WRITER _writer;
  
}; // class ParWriteSolution

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileWriter/ParWriteSolution.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileWriter_ParWriteSolution_hh

