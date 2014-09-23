// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileWriter_WriteSolution_hh
#define COOLFluiD_IO_CFmeshFileWriter_WriteSolution_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Table.hh"
#include "Framework/CFmeshFileWriter.hh"
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
/// @author Tiago Quintino
class CFmeshFileWriter_API WriteSolution : public CFmeshWriterCom {
public:

  typedef Framework::CFmeshFileWriter<Framework::CFmeshWriterSource> Writer;

  /// Constructor.
  explicit WriteSolution(const std::string& name);

  /// Destructor.
  ~WriteSolution();

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

  /// the file writer
  Writer m_writer;

  /// the actual data to write to the file
  std::auto_ptr<Framework::CFmeshWriterSource> m_data;

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;

}; // class WriteSolution

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileWriter_WriteSolution_hh

