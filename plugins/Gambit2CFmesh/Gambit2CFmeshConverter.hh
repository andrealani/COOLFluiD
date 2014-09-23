// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_Gambit2CFmesh_Gambit2CFmeshConverter_hh
#define COOLFluiD_IO_Gambit2CFmesh_Gambit2CFmeshConverter_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>




#include "Framework/MeshFormatConverter.hh"
#include "ElementTypeGambit.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace Gambit2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * Provides an abstract interface for the format converters of the mesh files
 * in input.
 *
 * @author Joao Pinto
 *
 */
class Gambit2CFmeshConverter : public Framework::MeshFormatConverter {
public:

  /**
   * Typedef of pair
   * 1. number = Gambit element type number
   * 2. number = number of nodes of the elemet
   */
  typedef std::pair<CFuint,CFuint> pairID;

  /**
   * This class stores information about the faces on the patches:
   */
  class FaceGambit {
  public:

    FaceGambit() :
      _faceNodes(0)
    {
    }

    CFuint getCellID() const
    {
      return _cellID;
    }

    CFuint& getCellID()
    {
      return _cellID;
    }

    void setCellID(const CFuint& cellID)
    {
      _cellID = cellID;
    }

    void setNbNodesInFace(const CFuint& nbNodesInFace)
    {
      cf_assert(nbNodesInFace > 0);
      _faceNodes.resize(nbNodesInFace);
    }

    CFuint getNbNodesInFace() const
    {
      return _faceNodes.size();
    }

    std::valarray<CFuint>& getFaceNodes()
    {
      return _faceNodes;
    }


  private:

    /// cell ID to which this face belongs
    CFuint _cellID;

    /// the IDs of th nodes in this face
    std::valarray<CFuint> _faceNodes;

  }; // end class FaceGambit

  /**
   * This class stores information about the topological patches:
   */
  class PatchGambit {
  public:

    PatchGambit() :
      _faceData(0)
    {
    }

    /**
     * Sets the Patch code
     */
    void setPatchCode(const CFuint& code)
    {
      _code = code;
    }

    /**
     * Sets the number of facesi the patch and
     *  resizes the storage for Face data.
     */
    void setNbFacesInPatch(const CFuint& nbFacesInPatch)
    {
      cf_assert(nbFacesInPatch > 0);
      _faceData.resize(nbFacesInPatch);
    }

    /**
     * Gets the number of Faces in the Patch
     */
    CFuint getNbFacesInPatch() const
    {
      return _faceData.size();
    }

    /**
     * Gets the Patch code
     */
    CFuint getPatchCode() const
    {
      return _code;
    }

    /**
     * Gets the Face data
     */
    std::vector<FaceGambit>& getFaceData()
    {
      return _faceData;
    }

  private:

    /// the code identifying this patch
    CFuint    _code;

    /// storage of the Faces Data
    std::vector<FaceGambit> _faceData;

  }; // end class PatchGambit

  /**
   * This class stores information about the Gambit SP (SP = super patches)
   */
  class SPdata {
  public:

    SPdata() :
      _patchIDs(0)
    {
    }

    void setSuperPatchName(const std::string& name)
    {
      _nameSP = name;
    }

    std::string getSuperPatchName() const
    {
      return _nameSP;
    }

    void setNbPatchesInSuperPatch(const CFuint& nbPatches)
    {
      _patchIDs.resize(nbPatches);
    }

    CFuint getNbPatchesInSuperPatch() const
    {
      return _patchIDs.size();
    }

    std::valarray<CFuint>& getPatchIDs()
    {
      return _patchIDs;
    }

  private:

    /// name that identifies the SuperPatch
    std::string _nameSP;

    /// list of the patch IDs composing
    /// this SuperPatch
    std::valarray<CFuint>  _patchIDs;

  }; // end class SPdata

public:

  /**
   * Constructor
   */
  Gambit2CFmeshConverter(const std::string& name);

  /**
   * Destructor
   */
  ~Gambit2CFmeshConverter();

  /**
   * Tries to check the file for conformity to the format.
   * Possibly not full proof.
   */
  void checkFormat(const boost::filesystem::path& filepath);

  /**
   * Writes the data read to the original format.
   * Useful for debugging purposes.
   */
  void convertBack(const boost::filesystem::path& filepath);

protected:

  /**
   * Adjust the node (state) numbering to make it stick to the
   * COOLFluiD convention.
   */
  void adjustToCFmeshNodeNumbering() {}

  /**
   * Gets the dimension.
   */
  CFuint getDimension() const {
    return _dimension;
  }

  /**
   * Helper function that throws an excetpion
   * when a error is found in the format of the file
   */
  void callGambitFileError(const std::string& msg,
       const CFuint& line,
       const boost::filesystem::path& filepath);

  /**
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getGambitWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);

  /**
   * Helper function that reads a line of cell/element IDs that belong
   * to a group in a Gambit neutral file and puts into a string incrementing
   * a supplied counter. Different from getGambitWordsFromLine to be able
   * to read meshes above 10 million cells/elements
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getGroupWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);

  /**
   * Helper function that reads the first line of connectivity data that belongs
   * to a cell/element in a Gambit neutral file and puts into a string incrementing
   * a supplied counter. Different from getGambitWordsFromLine to be able
   * to read meshes above 10 million cells/elements
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getFirstConnectivityWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);

  /**
   * Helper function that reads the other lines of connectivity data that belong
   * to a cell/element in a Gambit neutral file and puts into a string incrementing
   * a supplied counter. Different from getGambitWordsFromLine to be able
   * to read meshes above 10 million cells/elements
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getRestConnectivityWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);


  /**
   * Gets the target format.
   */
  std::string getTargetFormat() const
  {
    return "CFmesh";
  }

  /**
   * Gets the origin format.
   */
  std::string getOriginFormat() const
  {
    return "neu";
  }

  /**
   * Reads the all Dpl files and assembles the data in the converter.
   * It is always in a convert or convertBack, but only exectuted once.
   *
   * @throw Common::FilesystemException if a file cannot be open
   * @throw BadFormatException if a file is ill formated
   */
  void readFiles(const boost::filesystem::path& filepath);

private:

  /**
   * Reads parameters and stores them in atributes of class Gambit2CFmeshConverter
   */
  void readParameters(std::ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb);

  /**
   * Reads node coordinates and stores them in atributes of class Gambit2CFmeshConverter
   */
  void readNodes(std::ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb);

  /**
   * Reads elements of each group and stores them
   */
  void readGroups(std::ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb);

  /**
   * Reads element connectivity and stores it in atributes of class Gambit2CFmeshConverter
   */
  void readConnectivity(std::ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb);

  /**
   * Reads boundary conditions and stores them in atributes of class Gambit2CFmeshConverter
   */
  void readBC(std::ifstream& fin, const boost::filesystem::path& filepath, CFuint& lineNb, const CFuint& patchCode);

  /**
   * Reads the Gambit neutral file
   *
   * @pre the extension of the mesh file is ".neu"
   * @throw Common::FilesystemException if the file cannot be open
  */
  void readGambitFile(const boost::filesystem::path& filepath);

 /**
   * Reads the SuperPatch-file to know about the InnerCells
   *
   * @pre the extension of the mesh file is ".SP"
   * @throw Common::FilesystemException if the file cannot be open
   * @throw BadFormatException if the file is ill formated
   */
  void readSPInnerCells(const boost::filesystem::path& filepath);

  /**
   * Write in the COOLFluiD format the element list for a FEM mesh
   */
  void writeContinuousElements(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the state list for a FEM mesh
   */
  void writeContinuousStates(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the Topological Region Set data
   * considering to have FEM
   */
  void writeContinuousTrsData(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the mesh data considering
   * to have cell center FVM
   */
  void writeDiscontinuousElements(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the Topological Region Set data
   * considering to have FVM
   */
  void writeDiscontinuousTrsData(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the state list for a
   * cell centered FVM mesh
   */
  void writeDiscontinuousStates(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the node list
   */
  void writeNodes(std::ofstream& fout);

  /**
   * Get the number of super patches
   */
  CFuint getNbSuperPatches() const
  {
    return _superPatch.size();
  }

  /**
   * Get the number of element types
   */
  CFuint getNbElementTypes() const
  {
    return _elementType.size();
  }

  /**
   * Get the number of nodes per element type
   */
  CFuint getNbNodesPerElem(const CFuint typeID)
  {
    return _nodesPerElemTypeTable[typeID-1];
  }

  /**
  * Creates map connecting elements order to element type order
  * 1. column - element type order
  * 2. column - line in particular element type table
  */
  void createOrderChangeMap()
  {
    cf_assert(_nbCells > 0);
    _orderChangeMap = new Common::Table<CFuint>(_nbCells,2,std::numeric_limits<CFuint>::max());
  }

  /**
  * Gets the elements order map
  */
  Common::Table<CFuint>& getOrderChangeMap()
  {
    cf_assert(_orderChangeMap != CFNULL);
    return *_orderChangeMap;
  }

private:

  CFuint                          _dimension;

  CFuint                          _nbUpdatableNodes;

  CFuint                          _nbCells;

  CFuint                          _nbFaces;

  CFuint                          _nbPatches;

  CFuint                          _nbSuperPatches;

  /// storage of data about the ElementTypes
  /// including connectivity
  std::vector<ElementTypeGambit>  _elementType;

  /// storage of data about the patches
  std::vector<PatchGambit>        _patch;

  /// storage of data about the SuperPatches
  std::vector<SPdata>            _superPatch;

  std::valarray<CFuint>         _update;

  std::vector<CFreal>             _coordinate;

  Common::Table<CFreal>*           _variables;

  bool                                  _isWithSolution;

  bool                                  _isFileRead;

  CFuint &                               _nbUpdatableStates;

  /// Storage of the number of nodes per Gambit typeID
  std::vector<CFuint>            _nodesPerElemTypeTable;

  CFuint                                  _inFieldSP;

  /// variables added by Joao Pinto 06-12-2005:

  CFuint                           _nbGroups;

  CFuint                           _nbVelComponents;

  bool                          _isHighOrder;

  /// Map connecting elements order to element type order
  Common::Table<CFuint>*                                       _orderChangeMap;

  /// Tables of edge or face definitions
  std::vector< Common::Table<CFuint>* >               _faceDefinionTable;

  /// Table structure definitions
  std::vector< std::valarray<CFuint> >            _tablePattern;

  /// Map from Gambit element types to Face tables
  std::vector<pairID>                                         _gambitElemTypeID;

  /// Map from Gambit element types to Face tables
  std::map< pairID, Common::Table<CFuint>* >      _mapFaces;

  /// Map from Gambit order of element nodes to CF order of element nodes for Hexa
  std::map< CFuint, CFuint>                              _mapNodeOrderHexa;

  /// Map from Gambit order of element nodes to CF order of element nodes for Pyram
  std::map< CFuint, CFuint>                              _mapNodeOrderPyram;
  
  /// global connectivity table 
  Common::ConnectivityTable<CFuint> _tableConnectivity;
  
  // list of elements for each group
  std::vector< std::vector<CFuint> > _groupCellIDs;

  // names of the groups
  std::vector< std::string > _groupNames;

}; // end class Gambit2CFmeshConverter

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace Gambit2CFmesh

  } // end of namespace IO

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Gambit2CFmesh_Gambit2CFmeshConverter_hh
