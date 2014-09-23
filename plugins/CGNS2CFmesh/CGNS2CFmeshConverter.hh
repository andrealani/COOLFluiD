// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS2CFmesh_CGNS2CFmeshConverter_hh
#define COOLFluiD_CGNS2CFmesh_CGNS2CFmeshConverter_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Framework/MeshFormatConverter.hh"
#include "CGNS2CFmesh/ElementTypeCGNS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CGNS { class CGNSData; }

    namespace CGNS2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * Converts meshes from the CGNs format to the CFmesh format
 *
 * @author Tiago Quintino
 *
 */
class CGNS2CFmeshConverter : public Framework::MeshFormatConverter {

public: // nested classes

  /**
   * This class stores information about the faces on the patches:
   */
  class FaceCGNS {
  public: // funtctions

    FaceCGNS() :
      m_faceNodes(0)
    {
    }

    CFuint getCellID() const
    {
      return m_cellID;
    }

    CFuint& getCellID()
    {
      return m_cellID;
    }

    void setCellID(const CFuint& cellID)
    {
      cf_assert(cellID > 0);
      m_cellID = cellID;
    }

    void setNbNodesInFace(const CFuint& nbNodesInFace)
    {
      cf_assert(nbNodesInFace > 0);
      m_faceNodes.resize(nbNodesInFace);
    }

    CFuint getNbNodesInFace() const
    {
      return m_faceNodes.size();
    }

    std::valarray<CFuint>& getFaceNodes()
    {
      return m_faceNodes;
    }

  private: // data

    /// cell ID to which this face belongs
    CFuint m_cellID;

    /// the IDs of th nodes in this face
    std::valarray<CFuint> m_faceNodes;

  }; // end class FaceCGNS

  /**
   * This class stores information about the topological patches:
   */
  class PatchCGNS {
  public: // funtctions

    PatchCGNS() :
      m_faceData(0)
    {
    }

    /**
     * Sets the Patch code
     */
    void setPatchCode(const CFuint& code)
    {
      m_code = code;
    }

    /**
     * Sets the number of facesi the patch and
     *  resizes the storage for Face data.
     */
    void setNbFacesInPatch(const CFuint& nbFacesInPatch)
    {
      cf_assert(nbFacesInPatch > 0);
      m_faceData.resize(nbFacesInPatch);
    }

    /**
     * Gets the number of Faces in the Patch
     */
    CFuint getNbFacesInPatch() const
    {
      return m_faceData.size();
    }

    /**
     * Gets the Patch code
     */
    CFuint getPatchCode() const
    {
      return m_code;
    }

    /**
     * Gets the Face data
     */
    std::vector<FaceCGNS>& getFaceData()
    {
      return m_faceData;
    }

  private: // data

    /// the code identifying this patch
    CFuint    m_code;

    /// storage of the Faces Data
    std::vector<FaceCGNS> m_faceData;

  }; // end class PatchCGNS

  /**
   * This class stores information about the CGNS file.SP (SP = super patches)
   */
  class SPdata {
  public: // funtctions

    SPdata() : m_patchIDs(0)
    {
    }

    void setSuperPatchName(const std::string& name)
    {
      m_nameSP = name;
    }

    std::string getSuperPatchName() const
    {
      return m_nameSP;
    }

    void setNbPatchesInSuperPatch(const CFuint& nbPatches)
    {
      m_patchIDs.resize(nbPatches);
    }

    CFuint getNbPatchesInSuperPatch() const
    {
      return m_patchIDs.size();
    }

    std::valarray<CFuint>& getPatchIDs()
    {
      return m_patchIDs;
    }

  private: // data

    /// name that identifies the SuperPatch
    std::string m_nameSP;

    /// list of the patch IDs composing
    /// this SuperPatch
    std::valarray<CFuint>  m_patchIDs;

  }; // end class SPdata

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CGNS2CFmeshConverter(const std::string& name);

  /**
   * Destructor
   */
  ~CGNS2CFmeshConverter();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

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

protected: // functions

  /**
   * Gets the dimension.
   */
  CFuint getDimension() const {
    return m_dimension;
  }

  /**
   * Gets the number of variables
   */
  CFuint getNbVariables() const
  {
    return m_nbvariables;
  }

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
    return "cgns";
  }

protected: // helper functions

  /**
   * Reads the CGNS file
   *
   * @pre the extension of the mesh file is ".cgns"
   * @param filepath the full file path
   * @throw Common::FilesystemException if the file cannot be open
  */
  virtual void readFiles(const boost::filesystem::path& filepath);

private: // helper functions

  /**
   * Writes the THOR file when converting back to THOR
   * Useful for debugging purposes.
   * @see convertBack()
   */
  void writeTHOR(const boost::filesystem::path& filepath);

  /**
   * Writes the SP file when converting back to THOR
   * Useful for debugging purposes.
   * @see convertBack()
   */
  void writeSP(const boost::filesystem::path& filepath);

  /**
   * Adjust the node(= state) numbering to make it stick to the
   * COOLFluiD convention.
   */
  void adjustToCFmeshNodeNumbering();

  /**
   * Adjust the node(= state) numbering to make it stick to the
   * CGNS convention.
   */
  void adjustToCGNSNodeNumbering();

  /**
   * Adjust the node(= state) numbering to a certain offset
   * @see adjustToCFmeshNodeNumbering
   * @see adjustToCGNSNodeNumbering
   */
  void offsetNumbering(const CFint& offset);

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
    return m_superPatch.size();
  }

  /**
   * Get the number of element types
   */
  CFuint getNbElementTypes() const
  {
    return m_elementType.size();
  }

  /// @return the number of states in the type according to the order of the solution
  CFuint getNbStatesInType(const CFuint& typeID);

private: // data

  /// the storage for the CGNSData
  COOLFluiD::CGNS::CGNSData * m_cgns;

  /// String for configuring the order of the solution space
  std::string m_solOrderStr;

  /// Order of the solution space
  CFPolyOrder::Type m_solOrder;

  /// dimension of the mesh
  CFuint  m_dimension;

  /// nb of variables in cgns file
  CFuint  m_nbvariables;

  /// nb of cells in cgns file
  CFuint m_nbCells;

  /// nb of updatable nodes in cgns file
  CFuint m_nbUpdatableNodes;

  /// nb of boundary faces in cgns file
  CFuint m_nbFaces;

  /// nb of patches in cgns file
  CFuint  m_nbPatches;

  /// storage of data about the ElementTypes
  /// including connectivity
  std::vector<ElementTypeCGNS>     m_elementType;

  /// storage of data about the patches
  std::vector<PatchCGNS>           m_patch;

  /// storage of data about the SuperPatches
  std::vector<SPdata>              m_superPatch;

  std::valarray<CFuint>           m_update;

  /// connectivity cell to pair of node and state for the discontinuous CFPolyOrder::ORDER1 case
  Common::Table< std::pair<CFuint,CFuint> >* m_scon;

  bool                           m_isWithSolution;

  bool                           m_isFileRead;

  bool                           m_offsetCGNS;

}; // end class CGNS2CFmeshConverter

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNS2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS2CFmesh_CGNS2CFmeshConverter_hh
