// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileWriter_ParCFmeshFileWriter_hh
#define COOLFluiD_CFmeshFileWriter_ParCFmeshFileWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Framework/CFmeshWriterSource.hh"
#include "Common/SafePtr.hh"
#include "CFmeshFileWriter/CFmeshFileWriter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CFmesh format writer.
/// @author Andrea Lani
class CFmeshFileWriter_API ParCFmeshFileWriter : 
	public Config::ConfigObject {

public:

  /// Constructor.
  ParCFmeshFileWriter();

  /// Destructor.
  virtual ~ParCFmeshFileWriter();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Opens and starts to write to the given file.
  /// @throw Common::FilesystemException
  virtual void writeToFile(const boost::filesystem::path& filepath);

  /// Sets the pointer to the stored data
  void setWriteData(const Common::SafePtr<Framework::CFmeshWriterSource>& data)
  {
    _writeData = data;
  }
  
  /// Sets up private data
  void setup() {}
  
  /// Gets the objedct where data is stored
  Framework::CFmeshWriterSource& getWriteData()
  {
    cf_assert(_writeData.isNotNull());
    return * _writeData;
  }

  /// Get the file extension
  const std::string getWriterFileExtension() const
  {
    return std::string(".CFmesh");
  }

  /// Releases all temporary memory created while writing the file
  void releaseTemporaryWriteMemory()
  {
    getWriteData().releaseMemory();
  }

protected: // methods

  /// Writes to the given file.
  /// @throw Common::FilesystemException
  void writeToFileStream(const boost::filesystem::path& filepath,
  		 std::ofstream *const fout);

  /// Get the name of the reader
  const std::string getWriterName() const
  {
    return "ParCFmeshFileWriter";
  }

private: // helper functions

  /// Writes the extra variables info
  void writeExtraVarsInfo(std::ofstream *const fout);
  
  /// Writes the extra variables that are not node or state related
  void writeExtraVars(std::ofstream *const fout);

  /// Writes the version info in the CFmesh file
  void writeVersionStamp(std::ofstream *const fout);

  /// Writes the number of nodes, states, elements
  void writeGlobalCounts(std::ofstream *const fout);

  /// Writes the elements
  void writeElements(std::ofstream *const fout);

  /// Writes the element list
  void writeElementList(std::ofstream *const fout);

  /// Writes the list of nodes
  void writeNodeList(std::ofstream *const fout);

  /// Writes the list of state tensors and initialize the dofs
  void writeStateList(std::ofstream *const fout);

  /// Writes the all the data relative to all TRSs
  void writeTrsData(std::ofstream *const fout);

  /// Writes the all the geometric entities for the given TRS
  void writeGeoList(CFuint iTRS, std::ofstream *const fout);

protected: // data

  /// communicator
  MPI_Comm _comm;

  /// rank of this processor
  CFuint _myRank;

  /// number of processors
  CFuint _nbProc;

  /// I/O rank
  CFuint _ioRank;

  /// flag telling if the file is a new one
  CFuint _isNewFile;

  /// acquaintance of the data present in the CFmesh file
  Common::SafePtr<Framework::CFmeshWriterSource> _writeData;

  /// set keeping track of the files already created
  std::set<boost::filesystem::path> _fileList;

  /// set keeping track of the files already created
  std::map<boost::filesystem::path, long> _mapFileToStartNodeList;

}; // class ParCFmeshFileWriter

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileWriter_ParCFmeshFileWriter_hh

