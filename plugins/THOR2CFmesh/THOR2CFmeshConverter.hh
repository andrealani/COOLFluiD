// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_THOR2CFmesh_THOR2CFmeshConverter_hh
#define COOLFluiD_THOR2CFmesh_THOR2CFmeshConverter_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Framework/MeshFormatConverter.hh"
#include "THOR2CFmesh/ElementTypeTHOR.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace THOR2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * Provides an abstract interface for the format converters of the mesh files
 * in input.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class THOR2CFmeshConverter : public Framework::MeshFormatConverter {

public: // nested classes

  /**
   * This class stores information about the faces on the patches:
   */
  class FaceTHOR {
  public: // funtctions

    FaceTHOR() :
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

  }; // end class FaceTHOR

  /**
   * This class stores information about the topological patches:
   */
  class PatchTHOR {
  public: // funtctions

    PatchTHOR() :
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
    std::vector<FaceTHOR>& getFaceData()
    {
      return m_faceData;
    }

  private: // data

    /// the code identifying this patch
    CFuint    m_code;

    /// storage of the Faces Data
    std::vector<FaceTHOR> m_faceData;

  }; // end class PatchTHOR

  /**
   * This class stores information about the THOR file.SP (SP = super patches)
   */
  class SPdata {
  public: // funtctions

    SPdata() :
      m_patchIDs(0)
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
  THOR2CFmeshConverter(const std::string& name);

  /**
   * Destructor
   */
  ~THOR2CFmeshConverter();

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
    return m_nbLamVariables + m_nbTurbVariables;
  }

  /**
   * Helper function that throws an excetpion
   * when a error is found in the format of the file
   */
  void callTHORFileError(const std::string& msg,
       const CFuint& line,
       const boost::filesystem::path& filepath);

  /**
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getTHORWordsFromLine(std::ifstream& fin,
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
    return "thor";
  }

  /**
   * Reads the all Dpl files and assembles the data in the converter.
   * It is always in a convert or convertBack, but only exectuted once.
   *
   * @throw Common::FilesystemException if a file cannot be open
   * @throw BadFormatException if a file is ill formated
   */
  void readFiles(const boost::filesystem::path& filepath);

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
   * Reads the THOR file
   *
   * @pre the extension of the mesh file is ".thor"
   * @param filepath the full file path
   * @throw Common::FilesystemException if the file cannot be open
  */
  void readTHOR(const boost::filesystem::path& filepath);

  /**
   * Reads the SuperPatch-file
   *
   * @pre the extension of the mesh file is ".SP"
   * @param filepath the full file path
   * @throw Common::FilesystemException if the file cannot be open
   * @throw BadFormatException if the file is ill formated
   */
  void readSP(const boost::filesystem::path& filepath);

  /**
   * Adjust the node(= state) numbering to make it stick to the
   * COOLFluiD convention.
   */
  void adjustToCFmeshNodeNumbering();

  /**
   * Adjust the node(= state) numbering to make it stick to the
   * THOR convention.
   */
  void adjustToTHORNodeNumbering();

  /**
   * Adjust the node(= state) numbering to a certain offset
   * @see adjustToCFmeshNodeNumbering
   * @see adjustToTHORNodeNumbering
   */
  void offsetNumbering(const CFint& offset);

  /**
   * Checks the node numberings for a certain type of elements in the
   * lodaded data from teh thor files.
   */
  void checkNodeNumberings(const CFuint& dim,
                           ElementTypeTHOR& elemType,
                           const Common::Table<CFreal> *const nodes) const;

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

private: // helper functions

  /// @return the number of states in the type according to the order of the solution
  CFuint getNbStatesInType(const CFuint& typeID);

private: // data

  /// String for configuring the order of the solution space
  std::string m_solOrderStr;

  /// Order of the solution space
  CFPolyOrder::Type m_solOrder;

  /// dimension of the mesh
  CFuint  m_dimension;

  /// nb of laminar variables in thor
  CFuint  m_nbLamVariables;

  /// nb of turbulent variables in thor
  CFuint  m_nbTurbVariables;

  /// nb of cells in thor
  CFuint m_nbCells;

  /// nb of updatable nodes in thor
  CFuint m_nbUpdatableNodes;

  /// nb of faces in thor
  CFuint m_nbFaces;

  /// nb of patches in thor
  CFuint  m_nbPatches;

  /// storage of data about the ElementTypes
  /// including connectivity
  std::vector<ElementTypeTHOR>     m_elementType;

  /// storage of data about the patches
  std::vector<PatchTHOR>           m_patch;

  /// storage of data about the SuperPatches
  std::vector<SPdata>              m_superPatch;

  std::valarray<CFuint>           m_update;

  Common::Table<CFreal>*            m_coordinate;

  Common::Table<CFreal>*            m_variables;

  /// connectivity cell to pair of node and state for the discontinuous CFPolyOrder::ORDER1 case
  Common::Table< std::pair<CFuint,CFuint> >* m_scon;

  bool                           m_isWithSolution;

  bool                           m_isFileRead;

  bool                           m_offsetTHOR;

  CFuint &                        m_nbUpdatableStates;

}; // end class THOR2CFmeshConverter

//////////////////////////////////////////////////////////////////////////////

    } // namespace THOR2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_THOR2CFmesh_THOR2CFmeshConverter_hh
