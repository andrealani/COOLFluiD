// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileWriter_ParCFmeshBinaryFileWriter_hh
#define COOLFluiD_CFmeshFileWriter_ParCFmeshBinaryFileWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Config/ConfigObject.hh"
#include "Framework/CFmeshWriterSource.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Common/MPI/MPIError.hh"
#include "CFmeshFileWriter/CFmeshFileWriter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CFmesh binary format writer.
/// @author Andrea Lani
class CFmeshFileWriter_API ParCFmeshBinaryFileWriter : public Config::ConfigObject {

public:

  /// Constructor.
  ParCFmeshBinaryFileWriter();

  /// Destructor.
  virtual ~ParCFmeshBinaryFileWriter();

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
  void setup();

  /// Gets the object where data is stored
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

protected: // helper functions
  
  /// Writes to the given file.
  /// @throw Common::FilesystemException
  void writeToFileStream(const boost::filesystem::path& filepath);
  
  /// Get the name of the reader
  const std::string getWriterName() const 
  {
    return "ParCFmeshBinaryFileWriter";
  }
  
  /// Writes the extra variables info
  void writeExtraVarsInfo(MPI_File* fh);
  
  /// Writes the extra variables that are not node or state related
  void writeExtraVars(MPI_File* fh);
  
  /// Writes the version info in the CFmesh filex
  void writeVersionStamp(MPI_File* fh);

  /// Writes the number of nodes, states, elements
  void writeGlobalCounts(MPI_File* fh);
  
  /// Writes the elements
  void writeElements(MPI_File* fh);
  
  /// Writes the element list
  void writeElementList(MPI_File* fh);
  
  /// Writes the list of nodes
  void writeNodeList(MPI_File* fh);
  
  /// Writes the list of state tensors and initialize the dofs
  void writeStateList(MPI_File* fh);

  /// Writes the all the data relative to all TRSs
  void writeTrsData(MPI_File* fh);

  /// Writes the all the geometric entities for the given TRS
  void writeGeoList(CFuint iTRS, MPI_File* fh);
  
  /// Writes the end of the file
  void writeEndFile(MPI_File* fh);
  
  /// Write a key and a value
  template <typename T> 
  void writeKeyValue(MPI_File* fh, std::string key, bool onlyKey=true, T value=T())
  {
    cf_assert(key.size() < 30);
    if (key != "\n") {
      std::string buf(30,' ');         // initialize buffer with fixed size
      buf.replace(0, key.size(), key); // replace substring with given key
      MPI_File_write(*fh, &buf[0], (int)buf.size(), MPI_CHAR, &_status);
    }
    else {
      MPI_File_write(*fh, &key[0], (int)key.size(), MPI_CHAR, &_status);
    }
    
    if (!onlyKey) {
      T count = value;
      MPI_File_write(*fh, &count, 1, Common::MPIStructDef::getMPIType(&count), &_status);
    }
  }
  
  /// Write a buffer using all cores
  template <typename T> 
  void writeAll(const std::string& name, MPI_File* fh, MPI_Offset offset, T* buf, const CFuint bufSize) 
  {
    using namespace COOLFluiD::Common; 
    
    PE::Group& wgroup = PE::getGroup(0);
    const CFuint maxSendSize = _maxBuffSize/sizeof(T);
    MPI_Offset maxOffset = offset + bufSize*sizeof(T);
    CFLog(VERBOSE, "maxOffset   = " << maxOffset << "\n");
    CFLog(VERBOSE, "maxSendSize = " << maxSendSize << "\n");
    
    CFuint bufID = 0;
    CFuint bufLeft = bufSize;
    CFuint totBufLeft = bufLeft;
    MPI_Offset bufOff = offset; 
    do {
      const CFuint nbBuf =  bufLeft/maxSendSize;
      const CFuint wSize =  bufLeft%maxSendSize;
      cf_assert(wSize < maxSendSize);
      int wBufSize = (bufOff < maxOffset) ? ((nbBuf > 0) ? (int)maxSendSize : (int)wSize) : 0;
      cf_assert(wBufSize <= maxSendSize);
      cf_assert(wBufSize <= bufSize);
      cf_assert(bufID < bufSize);
      
      CFLog(VERBOSE, _myRank << " in " << name << " writes buffer of size " 
	    << wBufSize << "/" << bufSize << " starting from " << bufOff << "\n");
      CFLog(VERBOSE, "ELEM bufID    = " << bufID    << "\n");
      CFLog(VERBOSE, "ELEM wBufSize = " << wBufSize << "\n");
      CFLog(VERBOSE, "ELEM bufOff   = " << bufOff   << "\n");
      
      MPIError::getInstance().
	check("MPI_File_write_at_all", "ParCFmeshBinaryFileWriter::writeAll()",  
	      MPI_File_write_at_all(*fh, bufOff, &buf[bufID], wBufSize, MPIStructDef::getMPIType(&buf[bufID]), &_status));
      
      bufID   = std::min(bufID + (CFuint)wBufSize, bufSize); 
      bufLeft = std::max(bufLeft - (CFuint)wBufSize, (CFuint)0);
      bufOff  = std::min(bufOff + (MPI_Offset)(wBufSize*sizeof(T)), maxOffset);
      
      CFLog(VERBOSE, _myRank << " in " << name << " written buffer of size " << wBufSize << "/" << bufSize << " ending in " << bufOff << "\n");
      
      totBufLeft = 0.; 
      MPI_Allreduce(&bufLeft, &totBufLeft, 1, MPIStructDef::getMPIType(&totBufLeft), MPI_SUM, wgroup.comm);
      
      CFLog(VERBOSE, _myRank << " in " << name << " total buffer left " << totBufLeft << "\n");
    } while (totBufLeft > 0);
  }
  
  /// Test a buffer to check if it can be read correctly
  template <typename T>
  void testBuff(const std::string& name, MPI_Offset offset, T* buf, const CFuint bufSize) 
  {
    using namespace COOLFluiD::Common; 
    
    CFLog(INFO, "%%% ParCFmeshBinaryFileWriter::testBuff() for " << name << " %%% START\n");
    
    // a new file testIO.file is created
    // data are first written then read back and compared with original buffer
    PE::Group& wgroup = PE::getGroup(0);
    MPI_File fhh;
    MPI_File* fh = &fhh;
    
    // force offset to 0
    offset = 0;
    
    // write the buffer
    MPI_File_open(wgroup.comm, "testIO.file", MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, fh); 
    writeAll(name, fh, offset, buf, bufSize);
    MPI_File_close(fh);
    
    // read the buffer
    std::vector<T> readbuf(bufSize);
    MPI_File_open(wgroup.comm, "testIO.file", MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, fh); 
    readAll(name, fh, offset, &readbuf[0], bufSize);
    MPI_File_close(fh);
    
    CFuint count = 0;
    for (CFuint i = 0; i < readbuf.size(); ++i) {
      if (buf[i] != readbuf[i]) {
	CFLog(INFO, "#"<< ++count << " => buf[" << i << "] = " << buf[i] << " != "  
	      << "readbuf[" << i << "] = " << readbuf[i] << "\n"); 
      }
    }
    
    if (count > 0) exit(1);
    CFLog(INFO, "%%% ParCFmeshBinaryFileWriter::testBuff() for " << name << " %%% END\n");
  }
  
  template <typename T>
  void readAll(const std::string& name, MPI_File* fh, MPI_Offset offset, T* buf, const CFuint bufSize)
  {
    using namespace COOLFluiD::Common;
    
    const CFuint maxSendSize = _maxBuffSize/sizeof(T);
    MPI_Offset maxOffset = offset + bufSize*sizeof(T);
    CFuint bufID = 0;
    CFuint bufLeft = bufSize;
    CFuint totBufLeft = bufLeft;
    MPI_Offset bufOff = offset;
    do {
      const CFuint nbBuf =  bufLeft/maxSendSize;
      const CFuint wSize =  bufLeft%maxSendSize;
      int wBufSize = (bufOff < maxOffset) ? ((nbBuf > 0) ? (int)maxSendSize : (int)wSize) : 0;
      CFLog(VERBOSE, _myRank << " in " << name << " reads buffer of size " << wBufSize << "/" << bufSize << " starting from " << bufOff << "\n");
      cf_assert(bufID < bufSize);
      
      MPIError::getInstance().check
	("MPI_File_read_at_all", "ParCFmeshBinaryFileReader::readAll()",
	 MPI_File_read_at_all(*fh, bufOff, &buf[bufID], wBufSize, MPIStructDef::getMPIType(&buf[bufID]), &_status));
      
      bufID   = std::min(bufID + (CFuint)wBufSize, bufSize);
      bufLeft = std::max(bufLeft - (CFuint)wBufSize, (CFuint)0);
      bufOff  = std::min(bufOff + (MPI_Offset)(wBufSize*sizeof(T)), maxOffset);
      CFLog(VERBOSE, _myRank << " in " << name << " read buffer of size " << wBufSize << "/" << bufSize << " ending in " << bufOff << "\n");
      
      totBufLeft = 0.;
      MPI_Allreduce(&bufLeft, &totBufLeft, 1, MPIStructDef::getMPIType(&totBufLeft), MPI_SUM, _comm);
      CFLog(VERBOSE, _myRank << " in " << name << " total buffer left " << totBufLeft << "\n");
      
    } while (totBufLeft > 0);
  }
     
protected: // data
  
  /// Class holding all offsets defining the parallel file structure
  class Offset {
  public:
    /// start/end of the element list
    std::pair<MPI_Offset, MPI_Offset> elems;
    
    /// start/end of the nodes list
    std::pair<MPI_Offset, MPI_Offset> nodes;
    
    /// start/end of the states list
    std::pair<MPI_Offset, MPI_Offset> states;
    
    /// start/end of the TRS element list
    std::vector<std::pair<MPI_Offset, MPI_Offset> > TRS;
  };
  
  /// communicator
  MPI_Comm _comm;
  
  /// file handler
  MPI_File _fh;
  
  /// file status
  MPI_Status _status;
    
  /// rank of this processor
  CFuint _myRank;
  

  /// number of processors
  CFuint _nbProc;
  
  /// I/O rank
  CFuint _ioRank;
  
  /// group ID
  CFuint _myGroupID;
  
  /// offsets holder
  Offset _offset;
    
  /// flag telling if the file is a new one
  bool _isNewFile;

  /// flag telling if this processor is a writer
  bool _isWriterRank;
  
  /// acquaintance of the data present in the CFmesh file
  Common::SafePtr<Framework::CFmeshWriterSource> _writeData;

  /// set keeping track of the files already created
  std::set<boost::filesystem::path> _fileList;
  
  /// set keeping track of the files already created
  std::map<boost::filesystem::path, long long> _mapFileToStartNodeList;
  
  /// number of writers (and MPI groups)
  CFuint _nbWriters;
  
  /// maximu size of the buffer to write with MPI I/O
  int _maxBuffSize;
  
}; // class ParCFmeshBinaryFileWriter

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileWriter_ParCFmeshBinaryFileWriter_hh

