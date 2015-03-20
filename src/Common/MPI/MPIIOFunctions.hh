// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_MPI_MPIIOFunctions_hh
#define COOLFluiD_Common_MPI_MPIIOFunctions_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>
#include <string>

#include "Common/PE.hh"
#include "Common/CFLog.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Common/MPI/MPIError.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// The following class defines an interface to perform I/O operations in 
/// parallel
/// @author Andrea Lani
class MPIIOFunctions {
public:
  
  /// Write a key and a value
  template <typename T> 
  static void writeKeyValue(MPI_File* fh, std::string key, bool onlyKey=true, T value=T())
  {
    cf_assert(key.size() < 30);
    MPI_Status status;
    
    if (key != "\n") {
      std::string buf(30,' ');         // initialize buffer with fixed size
      buf.replace(0, key.size(), key); // replace substring with given key
      MPI_File_write(*fh, &buf[0], (int)buf.size(), MPI_CHAR, &status);
    }
    else {
      MPI_File_write(*fh, &key[0], (int)key.size(), MPI_CHAR, &status);
    }
    
    if (!onlyKey) {
      T count = value;
      MPI_File_write(*fh, &count, 1, Common::MPIStructDef::getMPIType(&count), &status);
    }
  }
  
  /// Write a buffer using all cores
  template <typename T> 
  static void writeAll(const std::string& name, MPI_File* fh, 
		       MPI_Offset offset, T* buf, 
		       const CFuint bufSize, 
		       int maxBuffSize) 
  {
    using namespace COOLFluiD::Common; 
    
    PE::Group& wgroup = PE::getGroup(0);
    const CFuint maxSendSize = maxBuffSize/sizeof(T);
    MPI_Offset maxOffset = offset + bufSize*sizeof(T);
    CFLog(VERBOSE, "maxOffset   = " << maxOffset << "\n");
    CFLog(VERBOSE, "maxSendSize = " << maxSendSize << "\n");
    
    int myRank = PE::GetPE().GetRank();
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
      
      CFLog(VERBOSE, myRank << " in " << name << " writes buffer of size " 
	    << wBufSize << "/" << bufSize << " starting from " << bufOff << "\n");
      CFLog(VERBOSE, "ELEM bufID    = " << bufID    << "\n");
      CFLog(VERBOSE, "ELEM wBufSize = " << wBufSize << "\n");
      CFLog(VERBOSE, "ELEM bufOff   = " << bufOff   << "\n");

      MPI_Status status;
      MPIError::getInstance().
	check("MPI_File_write_at_all", "MPIIOFunctions::writeAll()",  
	      MPI_File_write_at_all(*fh, bufOff, &buf[bufID], wBufSize, MPIStructDef::getMPIType(&buf[bufID]), &status));
      
      bufID   = std::min(bufID + (CFuint)wBufSize, bufSize); 
      bufLeft = std::max(bufLeft - (CFuint)wBufSize, (CFuint)0);
      bufOff  = std::min(bufOff + (MPI_Offset)(wBufSize*sizeof(T)), maxOffset);
      
      CFLog(VERBOSE, myRank << " in " << name << " written buffer of size " << wBufSize << "/" << bufSize << " ending in " << bufOff << "\n");
      
      totBufLeft = 0.; 
      MPI_Allreduce(&bufLeft, &totBufLeft, 1, MPIStructDef::getMPIType(&totBufLeft), MPI_SUM, wgroup.comm);
      
      CFLog(VERBOSE, myRank << " in " << name << " total buffer left " << totBufLeft << "\n");
    } while (totBufLeft > 0);
  }
  
  /// Test a buffer to check if it can be read correctly
  template <typename T>
  static void testBuff(const std::string& name, MPI_Offset offset, T* buf, const CFuint bufSize) 
  {
    using namespace COOLFluiD::Common; 
    
    CFLog(VERBOSE, "MPIIOFunctions::testBuff() for " << name << " START\n");
    
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
    CFLog(VERBOSE, "MPIIOFunctions::testBuff() for " << name << " END\n");
  }
  
  /// Read a scalar
  template <typename T>
  static void readScalar(MPI_File* fh, T& value)
  {
    MPI_Status status;
    MPI_File_read_all(*fh, &value, 1, Common::MPIStructDef::getMPIType(&value), &status);
  }
  
  /// Read an array
  template <typename T>
  static void readArray(MPI_File* fh, T* array, CFuint ns)
  {
    MPI_Status status;
    MPI_File_read_all(*fh, &array[0], (int)ns, Common::MPIStructDef::getMPIType(&array[0]), &status);
  }
  
  /// Read a buffer using all cores
  template <typename T>
  static void readAll(const std::string& name, MPI_File* fh, MPI_Offset offset, T* buf, const CFuint bufSize, int maxBuffSize)
  {
    const CFuint maxSendSize = maxBuffSize/sizeof(T);
    MPI_Offset maxOffset = offset + bufSize*sizeof(T);
    CFuint bufID = 0;
    CFuint bufLeft = bufSize;
    CFuint totBufLeft = bufLeft;
    MPI_Offset bufOff = offset;
    MPI_Status status;
    MPI_Comm comm = PE::GetPE().GetCommunicator();
    int myRank = PE::GetPE().GetRank();
    
    do {
      const CFuint nbBuf =  bufLeft/maxSendSize;
      const CFuint wSize =  bufLeft%maxSendSize;
      int wBufSize = (bufOff < maxOffset) ? ((nbBuf > 0) ? (int)maxSendSize : (int)wSize) : 0;
      CFLog(VERBOSE, myRank << " in " << name << " reads buffer of size " << wBufSize << "/" << bufSize << " starting from " << bufOff << "\n");
      cf_assert(bufID < bufSize);
      
      MPIError::getInstance().check
	("MPI_File_read_at_all", "ParCFmeshBinaryFileReader::readAll()",
	 MPI_File_read_at_all(*fh, bufOff, &buf[bufID], wBufSize, MPIStructDef::getMPIType(&buf[bufID]), &status));
      
      bufID   = std::min(bufID + (CFuint)wBufSize, bufSize);
      bufLeft = std::max(bufLeft - (CFuint)wBufSize, (CFuint)0);
      bufOff  = std::min(bufOff + (MPI_Offset)(wBufSize*sizeof(T)), maxOffset);
      CFLog(VERBOSE, myRank << " in " << name << " read buffer of size " << wBufSize << "/" << bufSize << " ending in " << bufOff << "\n");
      
      totBufLeft = 0.;
      MPI_Allreduce(&bufLeft, &totBufLeft, 1, MPIStructDef::getMPIType(&totBufLeft), MPI_SUM, comm);
      CFLog(VERBOSE, myRank << " in " << name << " total buffer left " << totBufLeft << "\n");
      
    } while (totBufLeft > 0);
  }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_MPI_MPIIOFunctions_hh
