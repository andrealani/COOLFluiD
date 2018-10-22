// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ParFileWriter_hh
#define COOLFluiD_Framework_ParFileWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>
#include <set>

#include "Framework/FileWriter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides a basic interface for parallel file writers
/// @author Andrea Lani
class Framework_API ParFileWriter : public FileWriter {
public:
  
  /// Default constructor without arguments
  ParFileWriter();
  
  /// Default destructor
  virtual ~ParFileWriter();
  
  /// Gets the file extension to append to the file name
  virtual const std::string getWriterFileExtension() const = 0;
  
 protected:
  
  /// Get the name of the writer
  virtual const std::string getWriterName() const = 0;
  
  /// Set group of writers
  virtual void setWriterGroup();  
  
  /// Set the default writers by picking an arbitrary selection of _nbWriters writers
  void setDefaultWriters(std::vector<int>& writerRanks);
  
  /// Set _nbWritersPerNode writers per node
  void setNodeWriters(std::vector<int>& writerRanks);
  
 protected: //data
  
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
  std::vector<Offset> _offset;
  
  /// flag telling if the file is a new one
  bool _isNewFile;

  /// flag telling if this processor is a writer
  bool _isWriterRank;
  
  /// set keeping track of the files already created
  std::set<boost::filesystem::path> _fileList;
  
  /// set keeping track of the files already created
  std::map<boost::filesystem::path, long long> _mapFileToStartNodeList;
  
  /// number of writers
  CFuint _nbWriters;
  
  /// number of writers per node
  CFuint _nbWritersPerNode;
  
  /// maximum size of the buffer to write with MPI I/O
  int _maxBuffSize;
  
  /// flag telling to write the FIRST CFmesh w/o solution
  bool _firstWithoutSolution;
  
}; // end of class ParFileWriter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ParFileWriter_hh
