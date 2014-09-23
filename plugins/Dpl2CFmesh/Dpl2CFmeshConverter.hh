// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_Dpl2CFmesh_Dpl2CFmeshConverter_hh
#define COOLFluiD_IO_Dpl2CFmesh_Dpl2CFmeshConverter_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Framework/MeshFormatConverter.hh"
#include "ElementTypeDpl.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace Dpl2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * Provides an abstract interface for the format converters of the mesh files
 * in input.
 *
 * @author Thomas Wuilbaut
 *
 */
class Dpl2CFmeshConverter : public Framework::MeshFormatConverter {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * This class stores information about the faces on the patches:
   */
  class FaceDpl {
  public:

    FaceDpl() :
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

  }; // end class FaceDpl

  /**
   * This class stores information about the topological patches:
   */
  class PatchDpl {
  public:

    PatchDpl() :
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
   * Sets the current index
   */
  void setCurrentIndex(const CFuint& index)
  {
    _index = index;
  }

  /**
   * Sets the index
   */
  CFuint getCurrentIndex()
  {
    return _index;
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
    std::vector<FaceDpl>& getFaceData()
    {
      return _faceData;
    }

  private:

    /// the code identifying this patch
    CFuint    _code;

    /// storage of the Faces Data
    std::vector<FaceDpl> _faceData;

    /// Face Index
    CFuint _index;

  }; // end class PatchDpl

  /**
   * This class stores information about the Dpl file.SP (SP = super patches)
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
  Dpl2CFmeshConverter(const std::string& name);

  /**
   * Destructor
   */
  ~Dpl2CFmeshConverter();

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
  void callDplFileError(const std::string& msg,
       const CFuint& line,
       const boost::filesystem::path& filepath);

  /**
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getDplWordsFromLine(std::ifstream& fin,
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
    return "dpl";
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
   * Reads the Dpl file
   *
   * @pre the extension of the mesh file is ".thor"
   * @throw Common::FilesystemException if the file cannot be open
  */
  void readDplFile(const boost::filesystem::path& filepath);

  /**
   * Reads the SuperPatch-file
   *
   * @pre the extension of the mesh file is ".SP"
   * @throw Common::FilesystemException if the file cannot be open
   * @throw BadFormatException if the file is ill formated
   */
  void readSPFile(const boost::filesystem::path& filepath);

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
   * Write in the COOLFluiD format the Topological Region Set data
   * considering to have FVM
   */
  void writeTrsDataCellCenterFVM(std::ofstream& fout);

  CFuint getNbSuperPatches() const
  {
    return _superPatch.size();
  }

  CFuint getNbElementTypes() const
  {
    return _elementType.size();
  }

  CFuint getNbNodesPerElem(const CFuint typeID)
  {
    return _nodesPerElemTypeTable[typeID-1];
  }

private:

  CFuint                           _dimension;

  CFuint                          _nbCells;

  CFuint                          _nbUpdatableNodes;

  CFuint                          _nbFaces;

  CFuint                           _nbPatches;

  CFuint                           _nbSuperPatches;

  /// storage of data about the ElementTypes
  /// including connectivity
  std::vector<ElementTypeDpl>     _elementType;

  /// storage of data about the patches
  std::vector<PatchDpl>           _patch;

  /// storage of data about the SuperPatches
  std::vector<SPdata>              _superPatch;

  std::valarray<CFuint>          _update;

  Common::Table<CFreal>*            _coordinate;

  Common::Table<CFreal>*            _variables;

  bool                           _isWithSolution;

  bool                           _isFileRead;

  CFuint &                        _nbUpdatableStates;

  /// Storage of the number of nodes per Dpl typeID
  std::vector<CFuint>             _nodesPerElemTypeTable;

  bool                          _isHybrid;

}; // end class Dpl2CFmeshConverter

//////////////////////////////////////////////////////////////////////////////

    } // namespace Dpl2CFmesh

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Dpl2CFmesh_Dpl2CFmeshConverter_hh
