// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_ParCFmeshBinaryFileReader_hh
#define COOLFluiD_CFmeshFileReader_ParCFmeshBinaryFileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/MPI/MPIStructDef.hh"
#include "CFmeshFileReader/ParCFmeshFileReader.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a parallel binary CFmesh format reader.
/// @author Andrea Lani
class CFmeshFileReader_API ParCFmeshBinaryFileReader : public ParCFmeshFileReader {
  
public: // member functions

  /// Constructor.
  ParCFmeshBinaryFileReader();

  /// Destructor.
  virtual ~ParCFmeshBinaryFileReader();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  // static void defineConfigOptions(Config::OptionList& options);
  
  /// Read the given file. This is a template method
  /// @throw Common::FilesystemException
  virtual void readFromFile(const boost::filesystem::path& filepath);
  
  /// Sets up private data
  virtual void setup();
  
protected: // functions
  
  /// Implementation of the hook method for finishing the reading.
  /// Deallocates temporary data.
  virtual void finish();
  
  /// Read an entry in the .CFmesh file
  bool readString(MPI_File* fh);
  
  /// Get the name of the reader
  virtual const std::string getReaderName() const
  {
    return "ParCFmeshBinaryFileReader";
  }
  
  /// Read a string and trim it
  std::string readAndTrimString(MPI_File* fh)
  {
    std::string key(30, ' ');
    MPI_File_read_all(*fh, &key[0], key.size(), MPI_CHAR, &m_status);
    Common::StringOps::trim2(key); // remove leading (none) and trailing spaces
    CFLog(VERBOSE, "CFmesh key = <" << key << ">\n");
    return key;
  }
  
  /// Read a scalar
  template <typename T>
  void readScalar(MPI_File* fh, T& value)
  {MPI_File_read_all(*fh, &value, 1, Common::MPIStructDef::getMPIType(&value), &m_status);}
  
  /// Read an array
  template <typename T>
  void readArray(MPI_File* fh, T* array, CFuint ns)
  {MPI_File_read_all(*fh, &array[0], ns, Common::MPIStructDef::getMPIType(&array[0]), &m_status);}
  
 private: // typedefs
  
  typedef std::vector<Framework::ElementTypeData> ElemTypeArray;

  /// pointer to ReaderFunction
  typedef void (ParCFmeshBinaryFileReader::*ReaderFun)(MPI_File* fh);
  
  /// type that maps a string read in the File with a ReaderFunction
  typedef std::map<std::string,
                   ReaderFun,
                   std::less<std::string> > MapString2ReaderP;
  
private: // member functions
  
  /// Sets m_mapString2Reader, that maps a given string to a corresponding
  /// reader function
  virtual void setMapString2Readers();
  
  /// Reads the space dimension
  void readCFVersion(MPI_File* fh);

  /// Reads the space dimension
  void readSvnVersion(MPI_File* fh);

  /// Reads the space dimension
  void readCFmeshVersion(MPI_File* fh);

  /// Reads the space dimension
  void readDimension(MPI_File* fh);

  /// Reads the number of equations
  void readNbEquations(MPI_File* fh);

  /// Reads the number of nodes
  void readNbNodes(MPI_File* fh);

  /// Reads the nb of dofs state tensors and initialize with them the dofs
  void readNbStates(MPI_File* fh);

  /// Reads if the past states are present in the file
  void readStorePastStates(MPI_File* fh);

  /// Reads if the past nodes are present in the file
  void readStorePastNodes(MPI_File* fh);

  /// Reads if the past states are present in the file
  void readStoreInterStates(MPI_File* fh);

  /// Reads if the past nodes are present in the file
  void readStoreInterNodes(MPI_File* fh);

  /// Reads the nb of extra variables associated with nodes
  void readNbExtraNodalVars(MPI_File* fh);

  /// Reads the nb of extra variables associated with states
  void readNbExtraStateVars(MPI_File* fh);

  /// Reads the nb of extra variables
  void readNbExtraVars(MPI_File* fh);
  
  /// Reads the names of the extra variables associated with nodes
  void readExtraNodalVarNames(MPI_File* fh);

  /// Reads the names of the extra variables associated with states
  void readExtraStateVarNames(MPI_File* fh);

  /// Reads the names of the extra variables
  void readExtraVarNames(MPI_File* fh);

  /// Reads the strides of the extra variables associated with nodes
  void readExtraNodalVarStrides(MPI_File* fh);

  /// Reads the strides of the extra variables associated with states
  void readExtraStateVarStrides(MPI_File* fh);

  /// Reads the strides of the extra variables
  void readExtraVarStrides(MPI_File* fh);
  
  /// Reads the extra variables associated with states
  void readExtraVars(MPI_File* fh);

  /// Reads the nb of elements
  void readNbElements(MPI_File* fh);

  /// Reads the nb of element types
  void readNbElementTypes(MPI_File* fh);

  /// Reads the order of the polynomial representation of the geometry
  void readGeometricPolyOrder(MPI_File* fh);

  /// Reads the order of the polynomial representation of the solution
  void readSolutionPolyOrder(MPI_File* fh);

  /// Reads the Type of the polynomial representation of the geometry
  void readGeometricPolyType(MPI_File* fh);

  /// Reads the Type of the polynomial representation of the solution
  void readSolutionPolyType(MPI_File* fh);

  /// Reads the element types (CFGeoShape::Type)
  void readElementTypes(MPI_File* fh);

  /// Reads the nb of elements per type
  void readNbElementsPerType(MPI_File* fh);

  /// Reads the nb of nodes per type
  void readNbNodesPerType(MPI_File* fh);

  /// Reads the nb of dofs per type
  void readNbStatesPerType(MPI_File* fh);

  /// Reads the list of nodes
  void readNodeList(MPI_File* fh);

  /// Reads the list of state tensors and initialize the dofs
  void readStateList(MPI_File* fh);

  /// Reads the data concerning the elements
  void readElementList(MPI_File* fh);

  /// Reads the number of topological region sets and initialize the vector
  /// that will contain the all the topological region sets
  /// @pre the element list has been already read
  /// @pre Connection has been already been and set
  /// @pre some topological region sets have been already constructed
  ///      (INNER_CELLS and, in FVM, INNER_FACES)
  void readNbTRSs(MPI_File* fh);

  /// Reads the name of the current topological region sets
  void readTRSName(MPI_File* fh);

  /// Reads the number of topological regions in the current
  /// topological region  set
  void readNbTRs(MPI_File* fh);

  /// Reads the number of geometric entities in each topological
  /// region of the current topological region set
  void readNbGeomEnts(MPI_File* fh);

  /// Reads the type of geometric entity in the current topological
  /// region set
  /// @pre  the read string must be "Face", "Cell" (or "Edge" in the future)
  /// @post the read string is converted in the corresponding CFGeoEnt::Type
  ///       by the method m_getCFGeoEnt::Type()
  void readGeomType(MPI_File* fh);

  /// Reads all the lists of geometric entities, using them to build the
  /// corresponding topological region.
  /// Once that all the topological regions belonging to the current topological
  /// region set have been built, the topological region set itself is built.
  void readGeomEntList(MPI_File* fh);

  /// Reads the data for one Topological Region Set
  void readTRSData(Common::CFMultiMap<CFuint,CFuint>& mapNodeElemID,
		   CFuint iTRS,
		   MPI_File* fh);
  
  /// Reads the element list corresponding for the current rank
  void readElemListRank( Framework::PartitionerData& pdata, MPI_File* fh);

  /// Reads the number of groups in the mesh
  void readNbGroups(MPI_File* fh);

  /// Reads the name of the current group
  void readGroupName(MPI_File* fh);

  /// Reads the number of elements in the current group
  void readGroupElementNb(MPI_File* fh);

  /// Reads the elements belonging to the current group in the mesh
  void readGroupElementList(MPI_File* fh);
  
  /// Set the donor rank and local ID corresponding to the given ID
  void setRankLocalID(const std::vector<std::pair<CFuint, CFuint> >& ranges, 
		      const CFuint globalID, CFuint& donorRank, CFuint& donorLocalID) const
  {
    for (CFuint i = 0; i < ranges.size(); ++i) {
      if (globalID >= ranges[i].first && globalID <= ranges[i].second) {
	donorRank    = i; // donor rank 
	donorLocalID = globalID - ranges[i].first; // local ID in donor
	break;
      } 
    }
  }
  
  /// Set send counts and displacements
  void setSendCountDispl(const std::vector<CFuint>& ranks,
			 const CFuint stride,
			 std::vector<int>& sendCount,
			 std::vector<int>& sendDispl);
  
  /// Set receive displacements
  void setRecvDispl(const std::vector<int>& recvCount,
		    std::vector<int>& recvDispl);
  
  /// Get the local (nodes or states) data
  void getLocalData(const std::vector<CFreal>& buf,
		    const std::vector<std::pair<CFuint, CFuint> > ranges,
		    const std::vector<CFuint>& listIDs,
		    const CFuint nodeSize, 
		    std::vector<CFreal>& recvBuf);
  
  /// Create the nodes storage
  void createNodesAll(const std::vector<CFreal>& localNodesData, 
		      const std::vector<CFreal>& ghostNodesData, 
		      Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes);
  
  /// Create the states storage
  void createStatesAll(const std::vector<CFreal>& localStatesData, 
		       const std::vector<CFreal>& ghostStatesData, 
		       Framework::DataHandle<Framework::State*, Framework::GLOBAL> states);
  
  /// Set the ranges to distribute reading load between processes 
  void setReadingRanges(const CFuint total, 
			std::vector<std::pair<CFuint, CFuint> >& ranges, 
			std::vector<CFuint>& nbNodesPerProc);
  
 private: // data
  
  /// map each string with a corresponding pointer to member
  /// function
  MapString2ReaderP m_mapString2ReaderFun;
  
  /// file handler
  MPI_File m_fh;
  
  /// file status
  MPI_Status m_status;
  
}; // class ParCFmeshBinaryFileReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileReader_ParCFmeshBinaryFileReader_hh
