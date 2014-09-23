// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_ParCFmeshFileReader_hh
#define COOLFluiD_CFmeshFileReader_ParCFmeshFileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMultiMap.hh"
#include "Common/FilesystemException.hh"

#include "Framework/FileReader.hh"
#include "Framework/CFmeshReaderSource.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/PartitionerData.hh"
#include "Framework/ElementDataArray.hh"

#include "CFmeshFileReader/CFmeshFileReaderAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework 
  {
      class VarSetTransformer;
      class MeshPartitioner;
  }

  namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a parallel CFmesh format reader.
/// @author Andrea Lani
class CFmeshFileReader_API ParCFmeshFileReader : public Framework::FileReader,
						 public Config::ConfigObject {

public: // member functions

  /// Constructor.
  ParCFmeshFileReader();
  
  /// Destructor.
  virtual ~ParCFmeshFileReader();
  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the data from the supplied arguments.
  /// @param args the argument list to configure this object
  void configure ( Config::ConfigArgs& args );
  
  /// Sets up private data
  virtual void setup();
    
  /// Sets the pointer to the stored data
  void setReadData(const Common::SafePtr<Framework::CFmeshReaderSource>& data)
  {
    cf_assert(data.isNotNull());
    m_readData = data;
  }

  /// Gets the object where data is stored
  Framework::CFmeshReaderSource& getReadData()
  {
    cf_assert(m_readData.isNotNull());
    return *m_readData;
  }

  /// Get the file extension
  virtual const std::string getReaderFileExtension() const
  {
    static const std::string ext = ".CFmesh";
    return ext;
  }

  /// Releases all temporary memory created while reading the file
  void releaseTemporaryReadMemory()
  {
    getReadData().releaseMemory();
  }

  /// Set some additional state values
  /// @param useInitValues  array of flags telling if to use the given
  ///                       initial values to initialize each variable in the
  ///                       State's
  /// @param initValues     initial values to set in the State
  void setStateInitValues(const std::vector<bool>& useInitValues,
                          const std::vector<CFreal>& initValues,
                          const std::vector<CFuint>& initValuesIDs)
  {
    m_useInitValues = useInitValues;
    m_initValues    = initValues;
    m_initValuesIDs = initValuesIDs;
  }

protected: // functions

  /// Implementation of the hook method for finishing the reading.
  /// Deallocates temporary data.
  virtual void finish();

  /// Read an entry in the .CFmesh file
  virtual bool readString(std::ifstream& file);

  /// Get the file extension
  virtual const std::string getReaderTerminator() const
  {
    static const std::string terminator = "!END";
    return terminator;
  }

  /// Get the name of the reader
  virtual const std::string getReaderName() const
  {
    return "ParCFmeshFileReader";
  }

private: // typedefs

  typedef std::vector<Framework::ElementTypeData> ElemTypeArray;

  /// pointer to ReaderFunction
  typedef void (ParCFmeshFileReader::*ReaderFunction)(std::ifstream& fin);

  /// type that maps a string read in the File with a ReaderFunction
  typedef std::map<std::string,
                   ReaderFunction,
                   std::less<std::string> > MapString2Reader;
  
protected: // member functions

  /// configures the merging of the TRS's
  void configureTRSMerging();

  /// Get the element type ID
  /// @param nbElemPerType global number of elements per type
  /// @param elemID        global element ID
  CFuint getElementType(ElemTypeArray& elementType, CFuint elemID)
  {
    CFuint count = 0;
    for (CFuint iType = 0; iType < elementType.size(); ++iType) {
      count += elementType[iType].getNbElems();
      if (elemID < count) return iType;
    }
    cf_assert(false);
    return 0;
  }

  /// Get the element count for iType less than current type
  CFuint getNewGlobalElementID(std::vector<Framework::MeshElementType>& elementType,
    	       CFuint globalElemID)
  {
    CFuint count = 0;
    for (CFuint iType = 0; iType < elementType.size(); ++iType) {
      if (globalElemID < elementType[iType].elementCount + count) {
cf_assert((globalElemID - count) < elementType[iType].elementCount);
return (globalElemID - count);
      }
      count += elementType[iType].elementCount;
    }
    cf_assert(false);
    return 0;
  }
  
  /// Sets m_mapString2Reader, that maps a given string to a corresponding
  /// reader function
  virtual void setMapString2Readers();
  
  /// Reads the space dimension
  void readCFVersion(std::ifstream& fin);

  /// Reads the space dimension
  void readSvnVersion(std::ifstream& fin);

  /// Reads the space dimension
  void readCFmeshVersion(std::ifstream& fin);

  /// Reads the space dimension
  void readDimension(std::ifstream& fin);

  /// Reads the number of equations
  void readNbEquations(std::ifstream& fin);

  /// Reads the number of nodes
  void readNbNodes(std::ifstream& fin);

  /// Reads the nb of dofs state tensors and initialize with them the dofs
  void readNbStates(std::ifstream& fin);

  /// Reads if the past states are present in the file
  void readStorePastStates(std::ifstream& fin);

  /// Reads if the past nodes are present in the file
  void readStorePastNodes(std::ifstream& fin);

  /// Reads if the past states are present in the file
  void readStoreInterStates(std::ifstream& fin);

  /// Reads if the past nodes are present in the file
  void readStoreInterNodes(std::ifstream& fin);

  /// Reads the nb of extra variables associated with nodes
  void readNbExtraNodalVars(std::ifstream& fin);

  /// Reads the nb of extra variables associated with states
  void readNbExtraStateVars(std::ifstream& fin);

  /// Reads the nb of extra variables
  void readNbExtraVars(std::ifstream& fin);
  
  /// Reads the names of the extra variables associated with nodes
  void readExtraNodalVarNames(std::ifstream& fin);

  /// Reads the names of the extra variables associated with states
  void readExtraStateVarNames(std::ifstream& fin);

  /// Reads the names of the extra variables
  void readExtraVarNames(std::ifstream& fin);

  /// Reads the strides of the extra variables associated with nodes
  void readExtraNodalVarStrides(std::ifstream& fin);

  /// Reads the strides of the extra variables associated with states
  void readExtraStateVarStrides(std::ifstream& fin);

  /// Reads the strides of the extra variables
  void readExtraVarStrides(std::ifstream& fin);
  
  /// Reads the extra variables associated with states
  void readExtraVars(std::ifstream& fin);

  /// Reads the nb of elements
  void readNbElements(std::ifstream& fin);

  /// Reads the nb of element types
  void readNbElementTypes(std::ifstream& fin);

  /// Reads the order of the polynomial representation of the geometry
  void readGeometricPolyOrder(std::ifstream& fin);

  /// Reads the order of the polynomial representation of the solution
  void readSolutionPolyOrder(std::ifstream& fin);

  /// Reads the Type of the polynomial representation of the geometry
  void readGeometricPolyType(std::ifstream& fin);

  /// Reads the Type of the polynomial representation of the solution
  void readSolutionPolyType(std::ifstream& fin);

  /// Reads the element types (CFGeoShape::Type)
  void readElementTypes(std::ifstream& fin);

  /// Reads the nb of elements per type
  void readNbElementsPerType(std::ifstream& fin);

  /// Reads the nb of nodes per type
  void readNbNodesPerType(std::ifstream& fin);

  /// Reads the nb of dofs per type
  void readNbStatesPerType(std::ifstream& fin);

  /// Reads the list of nodes
  void readNodeList(std::ifstream& fin);

  /// Reads the list of state tensors and initialize the dofs
  void readStateList(std::ifstream& fin);

  /// Reads the data concerning the elements
  void readElementList(std::ifstream& fin);

  /// Reads the number of topological region sets and initialize the vector
  /// that will contain the all the topological region sets
  /// @pre the element list has been already read
  /// @pre Connection has been already been and set
  /// @pre some topological region sets have been already constructed
  ///      (INNER_CELLS and, in FVM, INNER_FACES)
  void readNbTRSs(std::ifstream& fin);

  /// Reads the name of the current topological region sets
  void readTRSName(std::ifstream& fin);

  /// Reads the number of topological regions in the current
  /// topological region  set
  void readNbTRs(std::ifstream& fin);

  /// Reads the number of geometric entities in each topological
  /// region of the current topological region set
  void readNbGeomEnts(std::ifstream& fin);

  /// Reads the type of geometric entity in the current topological
  /// region set
  /// @pre  the read string must be "Face", "Cell" (or "Edge" in the future)
  /// @post the read string is converted in the corresponding CFGeoEnt::Type
  ///       by the method m_getCFGeoEnt::Type()
  void readGeomType(std::ifstream& fin);

  /// Reads all the lists of geometric entities, using them to build the
  /// corresponding topological region.
  /// Once that all the topological regions belonging to the current topological
  /// region set have been built, the topological region set itself is built.
  void readGeomEntList(std::ifstream& fin);

  /// Reads the data for one Topological Region Set
  void readTRSData(Common::CFMultiMap<CFuint,CFuint>& mapNodeElemID,
    CFuint iTRS,
    std::ifstream& fin);

  /// Reads the element list corresponding for the current rank
  void readElemListRank( Framework::PartitionerData& pdata, std::ifstream& fin);

  /// Reads the number of groups in the mesh
  void readNbGroups(std::ifstream& fin);

  /// Reads the name of the current group
  void readGroupName(std::ifstream& fin);

  /// Reads the number of elements in the current group
  void readGroupElementNb(std::ifstream& fin);

  /// Reads the elements belonging to the current group in the mesh
  void readGroupElementList(std::ifstream& fin);

  /// Ineffective reading of the node list
  void emptyNodeListRead(std::ifstream& fin);

  /// Ineffective reading of the state list
  void emptyStateListRead(std::ifstream& fin);

 protected:
  
  /// Set the element distribution array
  void setElmDistArray(std::vector<Framework::PartitionerData::IndexT>& elmdist);
  
  /// Set the size arrays for element-node and element-state connectivities
  void setSizeElemVec(std::vector<Framework::PartitionerData::IndexT>& sizeElemNodeVec,
		      std::vector<Framework::PartitionerData::IndexT>& sizeElemStateVec);
  
  /// Move the element data to the right processes and build info about the
  /// overlap region
  void moveElementData(Framework::ElementDataArray<0>& localElem,
		       Framework::PartitionerData& pdata);

  /// Add an element
  void addElement(const Framework::PartitionerData& pdata,
		  const CFuint globalElemID,
		  const CFuint localElemID,
		  Framework::ElementDataArray<0>& tmpElem,
		  CFuint& elemSize);
  
  /// Sort the m_pdata.part array to prepare a total exchange.
  void sortPartVec(const std::vector<CFuint>& globalElemID,
		   std::vector<Framework::PartitionerData::IndexT>& part,
		   Common::CFMultiMap<CFint,CFint>& sortedPart);
  
  /// Set the array telling if nodes and states are local or not
  void setIsLocalNodeState(Framework::ElementDataArray<0>& localElem,
			   std::vector<bool>& isLocalNode,
			   std::vector<bool>& isLocalState,
			   std::vector<CFuint>& localNodeIDs,
			   std::vector<CFuint>& localStateIDs);
  
  /// Update the array telling if nodes and states are local or not
  void updateIsLocalNodeState(CFuint root,
			      Framework::ElementDataArray<0>& elem,
			      Framework::ElementDataArray<0>& overlapElem,
			      std::vector<bool>& isLocalNode,
			      std::vector<bool>& isLocalState,
			      std::vector<CFuint>& ghostNodeIDs,
			      std::vector<CFuint>& ghostStateIDs,
			      std::vector<CFuint>& newLocalNodeIDs,
			      std::vector<CFuint>& newLocalStateIDs,
			      std::vector<CFuint>& localNodeIDsToRemove,
			      std::vector<CFuint>& localStateIDsToRemove,
			      std::vector<bool>& isOverlap,
			      CFuint nOverlap);
  
  /// Set the mapping between the global and the local node (or state) ID
  void setMapGlobalToLocalID(const std::vector<CFuint>& localIDs,
			     const std::vector<CFuint>& ghostIDs,
			     Common::CFMap<CFuint,CFuint>& m);
  
  /// Set the mapping between the global nodeID and the local elementID
  void setMapNodeElemID(Framework::ElementDataArray<0>& localElem);
  
  /// Set the elements
  void setElements(Framework::ElementDataArray<0>& localElem);

  /// Check the validity of a degree of freedom
  void checkDofID(CFuint dofID, CFuint totCount)
  {
    if (dofID >= totCount) {
      CFLog(ERROR, "DofID: " << dofID << " >= totCount \n");
      abort();
    }
  }
  
  /// Check if the given value is an entry in the container
  template <typename ARRAY, typename T>
    bool hasEntry(const ARRAY& array, const T& value)
  {
    return binary_search(array.begin(), array.end(), value);
  }
  
  /// Flag telling if the elements have been built
  bool areElementsBuild() const
  {
    return ((m_localNodeIDs.size() > 0) || (m_ghostNodeIDs.size() > 0)) &&
      ((m_localStateIDs.size() > 0) || (m_ghostStateIDs.size() > 0)) ;
  }
  
 private: // data
  
  /// map each string with a corresponding pointer to member
  /// function
  MapString2Reader m_mapString2Reader;
  
 protected:
  
  /// acquaintance of the data present in the CFmesh file
  Common::SelfRegistPtr<Framework::MeshPartitioner> m_partitioner;
  
  /// communicator
  MPI_Comm m_comm;

  /// local element storage
  Framework::ElementDataArray<0>* m_local_elem;

  /// rank of this processor
  CFuint m_myRank;

  /// number of processors
  CFuint m_nbProc;

  /// starting position of node list
  long unsigned int m_startNodeList;

  /// starting position of state list
  long unsigned int m_startStateList;

  /// array of flags telling if to use the given initial
  /// values to initialize each variable in the State's
  std::vector<bool> m_useInitValues;

  /// array of initial values to set in the State's
  std::vector<CFreal> m_initValues;

  /// array of initial values IDs to set in the State's
  std::vector<CFuint> m_initValuesIDs;

  /// acquaintance of the data present in the CFmesh file
  Common::SafePtr<Framework::CFmeshReaderSource> m_readData;

  /// original number of equations in the CFmesh file
  CFuint m_originalNbEqs;

  /// global number of elements
  CFuint m_totNbElem;

  /// global number of element types
  CFuint m_totNbElemTypes;

  /// global number of nodes
  CFuint m_totNbNodes;

  /// global number of states
  CFuint m_totNbStates;

  /// array storing the global number of elements per processor
  std::vector<CFuint> m_nbElemPerProc;
  
  /// array storing the coloring data coming from the mesh partitoner
  std::vector<Framework::PartitionerData::IndexT> m_partitionerOutData;
  
  /// local node IDs
  std::vector<CFuint> m_localNodeIDs;

  /// local state IDs
  std::vector<CFuint> m_localStateIDs;

  /// ghost node IDs
  std::vector<CFuint> m_ghostNodeIDs;

  /// ghost state IDs
  std::vector<CFuint> m_ghostStateIDs;

  /// map global node ID to local node ID
  Common::CFMap<CFuint,CFuint> m_mapGlobToLocNodeID;

  /// map global state ID to local state ID
  Common::CFMap<CFuint,CFuint> m_mapGlobToLocStateID;

  /// map nodeID to local element ID
  Common::CFMultiMap<CFuint,CFuint> m_mapNodeElemID;

  /// local element IDs
  std::vector<CFuint> m_localElemIDs;

  ///PastStates and PastNodes present in CFmesh?
  bool m_hasPastNodes;
  bool m_hasPastStates;

   ///InterStates and InterNodes present in CFmesh?
  bool m_hasInterNodes;
  bool m_hasInterStates;

  /// name of the TRS currently being read
  /// used for merging the TRS's
  std::string m_curr_trs;

  /// reverse map from TRS to be merged to summed TRS
  std::map<std::string,std::string> m_rev_trs_merge;

  /// map from TRS to build to the index of that TRS
  std::map<std::string,CFuint> m_trs_idxmap;

  /// change in number of TRS's
  CFuint m_trs_reduction;

  /// merged TRS just added
  bool m_mergedtrs_just_added;
  /// current TRS is being merged
  bool m_mergedtrs;

  /// current number of TRs
  CFuint m_curr_nbtr;

  /// number of layers of overlap region
  CFuint m_nbOverLayers;

  /// partitioner name
  std::string m_partitionerName;

  /// config option for merging th TRS's
  std::vector<std::string> m_merge_trs;

  /// map from Group name to build to the index of that Group
  std::map<std::string,CFuint> m_groups_idxmap;

  /// name of the Group currently being read
  std::string  m_curr_group;

  /// Name of the vector transformer from input to update variables
  std::string m_inputToUpdateVecStr;

  /// Vector transformer from input to update variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVecTrans;

}; // class ParCFmeshFileReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileReader_ParCFmeshFileReader_hh
