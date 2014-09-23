// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Tecplot2CFmesh_Tecplot2CFmeshConverter_hh
#define COOLFluiD_Tecplot2CFmesh_Tecplot2CFmeshConverter_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Common/CFMap3D.hh"
#include "Framework/MeshFormatConverter.hh"
#include "Tecplot2CFmesh/ElementTypeTecplot.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Tecplot2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * Provides an abstract interface for the format converters of the mesh files
 * in input.
 *
 * @author Andrea Lani
 *
 */
class Tecplot2CFmeshConverter : public Framework::MeshFormatConverter {

public: // functions
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  Tecplot2CFmeshConverter(const std::string& name);

  /**
   * Destructor
   */
  ~Tecplot2CFmeshConverter();

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
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getTecplotWordsFromLine(std::ifstream& fin, std::string& line,
			       CFuint&  lineNb, std::vector<std::string>& words)
  {
    getline(fin,line);
    ++lineNb;
    words = Common::StringOps::getWords(line);
  }
  
  /**
   * Gets the target format.
   */
  std::string getTargetFormat() const {return "CFmesh";}

  /**
   * Gets the origin format.
   */
  std::string getOriginFormat() const {return "plt";}
  
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
   * Writes the Tecplot file when converting back to Tecplot
   * Useful for debugging purposes.
   * @see convertBack()
   */
  void writeTecplot(const boost::filesystem::path& filepath);
  
  /**
   * Reads the individual Tecplot file
   *
   * @param filepath the full file path
   * @throw Common::FilesystemException if the file cannot be open
   */
  void readTecplotFile(CFuint nbZones, 
		       const boost::filesystem::path& filepath, 
		       const std::string& extension, 
		       bool isBoundary);
  
  /**
   * Adjust the node(= state) numbering to a certain offset
   */
  void offsetNumbering(CFint offset);

  /// Adjust the node (state) numbering to make it stick to the
  /// COOLFluiD convention.
  void adjustToCFmeshNodeNumbering();
  
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
   * Get the number of element types
   */
  CFuint getNbElementTypes() const
  {
    return m_nbElemTypes;
  }
  
  /**
   * Get the number of nodes per cell
   */
  CFuint getNbNodesInCell(const std::string& etype) const;  
  
  /**
   * Get the space dimension
   */
  CFuint getDim(const std::string& etype) const;

private: // helper functions

  /// @return the number of states in the type according to the order of the solution
  CFuint getNbStatesInType(const CFuint& typeID);

  // get the value string corresponding to the given key
  void getValueString(const std::string& key, const std::string& in, std::string& out) const;
  
  // read the variables information
  void readVariables(std::ifstream& fin, std::string& line, 
		     CFuint&  lineNb, 
		     std::vector<std::string>& words, 
		     std::vector<std::string>& vars);
  
  // read element information (mesh/solution nodes, connectivity) from a zone
  void readZone(std::ifstream& fin, CFuint& countl, 
		std::vector<std::string>& vars,
		bool isBoundary);
  
  // renumber tthe TRS data
  template <CFuint DIM> 
  void renumberTRSData();
  
  // build the connectivity between boundary faces and neighboring cells
  void buildBFaceNeighbors();
  
  // skip the solution
  bool skipSolution() const {return (m_readVars.size() == 0);}
  
private: // data

  /// dimension of the mesh
  CFuint  m_dimension;
   
  /// storage of data about the ElementTypes
  /// including connectivity
  std::vector<ElementTypeTecplot*> m_elementType;
    
  /// data organized in blocks
  bool m_hasBlockFormat;
  
  /// surface TRS names
  std::vector<std::string>  m_surfaceTRS; 
  
  /// variables to read
  std::vector<std::string>  m_readVars; 
  
  /// number of element types
  CFuint  m_nbElemTypes;
    
}; // end class Tecplot2CFmeshConverter

//////////////////////////////////////////////////////////////////////////////

template <CFuint DIM>
void Tecplot2CFmeshConverter::renumberTRSData()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Environment;
  
  // assume 1 element type
  SafePtr<ElementTypeTecplot> elements = m_elementType[0];
  const CFuint nbElems = elements->getNbElems();  
  const CFuint nbNodesPerElem = elements->getNbNodesPerElem();  
  
  // you could put a loop around this to analyze blocks of data instead of the whole in one shot
  typedef NodeDim<DIM> NodeD;
  CFMultiMap<NodeD, CFuint> mapCoord2nodeIDs(elements->getNbNodes());
  
  cout << "Tecplot2CFmeshConverter::renumberTRSData(), DIM, NBVARS = " << m_dimension << ", " << elements->getNbVars() << endl;
  
  vector<bool> flagNode(elements->getNbNodes(), false);
  NodeD nodeDim(CFNULL);
  for (CFuint e = 0; e < nbElems; ++e) {
    for (CFuint iNode = 0; iNode < nbNodesPerElem; ++iNode) {
      assert(e < elements->getNbElems());
      assert(iNode < nbNodesPerElem);
      const CFuint nodeID = (*elements)(e, iNode);
      if (nodeID >= elements->getNbNodes()) {
	cout << nodeID << " >= " << elements->getNbNodes() << endl;
	assert(nodeID < elements->getNbNodes());
      }
      
      if (!flagNode[nodeID]) {
	elements->setNodeDim(nodeID, nodeDim);
	mapCoord2nodeIDs.insert(nodeDim, nodeID);
	flagNode[nodeID] = true;
      }
    }
  }
  mapCoord2nodeIDs.sortKeys();
  
  cout << "Tecplot2CFmeshConverter::renumberTRSData(), insertion finished\n";
  
  typedef typename CFMultiMap<NodeD, CFuint>::MapIterator MapIt;
  
  CFuint counter = 0;
  CFuint totalNbNodes = 0;
  const CFuint startTRS = m_nbElemTypes;
  const CFuint nbTRSs = m_elementType.size() - startTRS;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<ElementTypeTecplot> trs = m_elementType[startTRS +  iTRS];
    const CFuint nbTRSFaces = trs->getNbElems();
    const CFuint nbNodesPerFace = trs->getNbNodesPerElem();
    for (CFuint iFace = 0; iFace < nbTRSFaces; ++iFace) {
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
	totalNbNodes++;
 	const CFuint trsNodeID = (*trs)(iFace,iNode);
	trs->setNodeDim(trsNodeID, nodeDim);

	bool found = false;
	pair<MapIt, MapIt> nodeList = mapCoord2nodeIDs.find(nodeDim, found);
	if (!found) cout << "Tecplot2CFmeshConverter::renumberTRSData() => nodeDim " << nodeDim << " not found!\n";
	
	// find(nodeDim, found);
	//if (found) {
	CFuint nbMatch = 0;
	for (MapIt it = nodeList.first; it != nodeList.second; ++it) {
	  (*trs)(iFace,iNode) = it->second;
	  nbMatch++;
	}
	// only one match allowed
	if (nbMatch != 1) {
	  cout << "Tecplot2CFmeshConverter::renumberTRSData() => " << nbMatch << " matches found for node (" << nodeDim << ")\n";
	  if (nbMatch == 0) {
	    cout << "Tecplot2CFmeshConverter::renumberTRSData() => node (" << nodeDim << ") not found" << endl;
	    counter++;
	  }
	  assert(nbMatch == 1);
	}
      }
    }
  }
  
  cout << "Tecplot2CFmeshConverter::renumberTRSData() => " << counter << "/" << totalNbNodes << " NOT FOUND" << endl;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Tecplot2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Tecplot2CFmesh_Tecplot2CFmeshConverter_hh
