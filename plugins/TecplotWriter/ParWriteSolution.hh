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

 extern "C" {
   static void cmpAndTakeMaxAbs3(CFreal* invec, CFreal* inoutvec, int* len,
				 MPI_Datatype* datatype)
   {
     cf_assert(len != CFNULL);
     int size = *len;
     for (int i = 0; i < size; ++i) {
       inoutvec[i] = (fabs(invec[i]) > 0.) ? invec[i] : inoutvec[i];
     }
   }
 }
      
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
  virtual ~ParWriteSolution();
  
  /// Set up private data
  virtual void setup();
  
  /// Unsetup private data
  virtual void unsetup();
  
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
  
  /// typedef for pointer to member function
  typedef void (ParWriteSolution::*WriterFun)
    (const boost::filesystem::path&, const bool, const std::string, std::ofstream*);
  
  /// This class stores indexes and offsets useful for the writing TRS data 
  class TecplotTRSType {
  public:
    /// Constructor
    TecplotTRSType(const CFuint nbElementTypes, 
		  Common::SafePtr<Framework::TopologicalRegionSet> in);
    
    /// Destructor
    ~TecplotTRSType();
    
    /// Clean up mappings
    void cleanupMappings();
    
    /// pointer to the corresponding TRS
    Common::SafePtr<Framework::TopologicalRegionSet> trs;
    
    /// start/end of the header
    std::vector<std::vector<MPI_Offset> > headerOffset;
    
    /// start/end of the element list
    std::vector<std::pair<MPI_Offset, MPI_Offset> > elemsOffset;
    
    /// start/end of the nodes list
    std::vector<std::pair<MPI_Offset, MPI_Offset> > nodesOffset;
    
    /// old total number of nodes and elements in type
    std::vector<std::pair<CFuint, CFuint> > oldNbNodesElemsInType;
    
    /// total number of nodes in element type
    std::vector<CFuint> totalNbNodesInType;
    
    /// global node IDs per type
    std::vector<std::vector<CFuint> > nodesInType;
    
    /// mapping global nodeIDs to global nodeIDs by element type
    std::vector<Common::CFMap<CFuint, CFuint>*> mapNodeID2NodeIDByEType;
  };
  
  /// write unstructured data on the given file
  /// @param filepath  namem of the path to the file
  /// @param flag      flag telling if the file has to be created or overwritten
  /// @param fun       pointer to member function to be called for writing
  void writeData(const boost::filesystem::path& filepath, 
		 bool& flag, 
		 const std::string title, 
		 WriterFun fun);
  
  /// Writes the inner data to file
  /// @param filename  name of output file
  /// @param isNewFile flag to tell if the file is to be created or to be overwritten
  /// @param fout      pointer to the file  
  virtual void writeInnerData(const boost::filesystem::path& filepath,
			      const bool isNewFile,
			      const std::string title,
			      std::ofstream* fout);
  
  /// Writes the boundary data to file
  /// @param filename  name of output file
  /// @param isNewFile flag to tell if the file is to be created or to be overwritten
  /// @param fout      pointer to the file  
  virtual void writeBoundaryData(const boost::filesystem::path& filepath,
				 const bool isNewFile,
				 const std::string title,
				 std::ofstream* fout);
  
  /// Writes the TECPLOT header
  // void writeHeader(MPI_File* fh);
  
  /// Writes the TECPLOT header
  virtual void writeHeader(std::ofstream* fout, 
			   const std::string title,
			   Common::SafePtr<Framework::DataHandleOutput> dh);
  
  /// Writes the TECPLOT zone header
  virtual void writeZoneHeader(std::ofstream* fout, 
			       const CFuint iType,
			       const std::string& geoShape,
			       const CFuint nbNodesInType,
			       const CFuint nbElemsInType,
			       const std::string& geoType,
			       const std::string& end,
			       const bool isBoundary); 
  
  /// Writes the TECPLOT inner zone header
  virtual void writeInnerZoneHeader(std::ofstream* fout, 
				    const CFuint iType,
				    Framework::ElementTypeData& eType,
				    Common::SafePtr<Framework::TopologicalRegionSet> trs); 
  
  /// Write the Tecplot file in Binmary format
  /// @throw Common::FilesystemException
  virtual void writeToBinaryFile();
  
  /// Write the node list corresponding to the given element type
  virtual void writeNodeList(std::ofstream* fout, const CFuint iType, 
			     Common::SafePtr<Framework::TopologicalRegionSet> elements,
			     const bool isBoundary);
  
  /// Write the element list corresponding to the given element type
  void writeElementList(std::ofstream* fout,
			const CFuint iType,
			const CFuint nbNodesInType,
			const CFuint nbElementsInType,
			const CFuint nbLocalElements,
			const CFuint startElementID,
			const CFuint geoOrder,
			Common::SafePtr<Framework::TopologicalRegionSet> elements);
  
  /// Get the number of element nodes to write 
  CFuint getWriteNbNodesInType(const CFuint nbNodes,
			       const CFuint geoOrder,
			       const CFuint dim, 
			       const bool isCell);
  
  /// Write the connectivity for one element, after having translated the local 
  /// element connectivity (node ordering) from CFmesh to Tecplot format
  void writeElementConn(std::ofstream& file,
			CFuint* nodeIDs,
			const CFuint nbNodes,
			const CFuint geoOrder,
			const CFuint dim, 
			const bool isCell);
  
  /// Get the name of the writer
  const std::string getWriterName() const;
  
  /// Tell if the mesh has changed 
  bool hasChangedMesh(const TecplotTRSType& tt) const;    
  
  /// Tell if the mesh has changed for a certain mesh type
  /// @param iType  ID for the mesh type
  bool hasChangedMesh(const CFuint iType, const TecplotTRSType& tt) const;    
  
  /// Return the TRS ID corresponding to the given TRS name in the global storage
  int getGlobalTRSID(const std::string& name) const;
  
  /// Build the mapping between global nodeIDs and global nodeIDs by element type
  void buildNodeIDMapping();
  
  /// Cleanup the mapping between global nodeIDs and global nodeIDs by element type
  void cleanupNodeIDMapping();
  
  /// Store the mapping coresponding to the given TRS and iType (element type or TR)
  void storeMappings(Common::SafePtr<Framework::TopologicalRegionSet> trs,
		     const CFuint iType,
		     const CFuint nbCellsInType,
		     const CFuint nbNodesInType,
		     const CFuint startIdx, 
		     const CFuint endIdx,
		     std::vector<bool>& foundGlobalID);
  
protected:
 
  /// socket for Node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// socket for State Proxy
  Framework::DataSocketSink<Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;
  
  /// array of TecplotTRSType
  Common::CFMap<std::string, TecplotTRSType*> _mapTrsName2TecplotData;
  
  /// mapping global nodeIDs to local nodeIDs
  Common::CFMap<CFuint, CFuint> _mapGlobal2LocalNodeID;
  
  /// vector with the names of the nodal variables to output
  std::vector<std::string> m_nodalvars;
  
  /// vector with the names of the cell-centered variables to output
  std::vector<std::string> m_ccvars;
  
  /// maximum size to be used for formatting integers
  CFuint _intWordFormatSize;
  
  /// flag to tell if the boundary file is new
  bool _isNewBFile;
  
  /// flag that specifies to output cell-centered or nodal variables
  bool m_onlyNodal;
  
  /// File format to write in (ASCII or Binary)
  std::string _fileFormatStr;
    
}; // class ParWriteSolution

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_TecplotWriter_ParWriteSolution_hh

