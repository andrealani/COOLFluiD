// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_TecplotWriter_ParWriteSolution_hh
#define COOLFluiD_IO_TecplotWriter_ParWriteSolution_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ParFileWriter.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

#include "TecplotWriter/TecWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a parallel writer for TECPLOT files
/// @author Andrea Lani
class TecplotWriter_API ParWriteSolution : 
	public TecWriterCom,
	public Framework::ParFileWriter {
  
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit ParWriteSolution(const std::string& name);

  /// Destructor.
  virtual ~ParWriteSolution()
  {
  }

    /// Set up private data
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  virtual void execute();
  
  /// Gets the file extension to append to the file name
  const std::string getWriterFileExtension() const
  {
    return std::string(".plt");
  }
  
 protected:
  
  /// Writes to the given file.
  /// @throw Common::FilesystemException
  void writeToFileStream(const boost::filesystem::path& filepath);
  
  /// Writes the TECPLOT header
  void writeHeader(MPI_File* fh);
  
  /// Write the Tecplot file in Binmary format
  /// @throw Common::FilesystemException
  virtual void writeToBinaryFile();
  
  /// Write the boundary surface data
  virtual void writeBoundarySurface();
  
  /// Get the name of the writer
  const std::string getWriterName() const;

protected:
  
  /// socket for Node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// socket for State Proxy
  Framework::DataSocketSink<Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;
  
  //File format to write in (ASCII or Binary)
  std::string _fileFormatStr;

}; // class ParWriteSolution

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_TecplotWriter_ParWriteSolution_hh

