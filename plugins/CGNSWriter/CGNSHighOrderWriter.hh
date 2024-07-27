// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CGNSWriter_CGNSHighOrderWriter_hh
#define COOLFluiD_IO_CGNSWriter_CGNSHighOrderWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FileWriter.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

#include "CGNSWriter/CGWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

    namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to write the MeshData
/// solution to a CGNS format file for visualization using high-order elements.
/// The method is tailored for FR solution data as it automatically upgrades the
/// MeshData to high geomtric order (Q) if the polynomial order P > Q of the orignal MeshData.
/// @author Rayan Dhib

class CGNSWriter_API CGNSHighOrderWriter : public CGWriterCom,
                      public Framework::FileWriter {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit CGNSHighOrderWriter(const std::string& name);

  /// Destructor.
  ~CGNSHighOrderWriter()
  {
  }

    /// Set up private data
  void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

  /// Gets the file extension to append to the file name
  const std::string getWriterFileExtension() const
  {
    return std::string(".cgns");
  }

protected:

  /// Write the to the given file stream the MeshData.
  /// @throw Common::FilesystemException
  void writeToFile(const std::string& fileName);

  /// Write the boundary surface data
  void writeBoundarySurface();

  /// Get the name of the writer
  const std::string getWriterName() const;

  /// Opens the CGNS file for writing
  void openFile(const std::string& fileName);

  /// Closes the CGNS file
  void closeFile();

  /// Write CGNS Base
  void writeCGNSBase(const CFuint cellDim, const CFuint physDim, std::string BaseName, int& fileIndex, int& baseIndex);

  /// Write CGNS Zone
  void writeCGNSZone(cgsize_t size[3], std::string ZoneName, int& fileIndex, int& baseIndex, int& index_zone);

  /// Write nodal coordinates
  void writeCGNSNodeCoordinates(CFuint totNbNodes, CFuint dim, std::vector< RealVector > nodes,int& fileIndex, int& baseIndex, int& index_zone, int& index_coord);

  /// Write connectivity
  void writeCGNSConnectivity(const std::vector<std::vector<cgsize_t>>& Connectivity, ElementType_t elementType, CFuint nbrElems, CFuint nbrNodes, int start, int end, int& fileIndex, int& baseIndex, int& index_zone, int& index_section) ;

  /// Write CGNS solution node
  void writeCGNSSolutionNode(const std::string& solutionName, int& fileIndex, int& baseIndex, int& index_zone, int& index_sol);

  /// returns the local subcell node connectivity in a cell of given shape and order (in CGNS conventions)
  static std::vector<CFuint> getCGNSCellNodeConn(ElementType_t elementType);

  /// returns the CGNS cell shape
  static ElementType_t getCGNSCellShape(CFGeoShape::Type shape,CFuint geoOrder);

  /// returns the nodal coords in the reference domain depending on geoOrder in CGNS conventions
  std::vector< RealVector > getOutputPntsMappedCoords(ElementType_t type);

  /// returns the nodal coords of the additional nodes in the physical domain
  std::vector< RealVector > calculateNewNodePositionsForElement(std::vector< RealVector >  outputPntsMappedCoords, CFuint nbrNodes, std::vector<COOLFluiD::Framework::Node*>*& cellNodes, std::vector< RealVector > geoShapeFuncs, CFuint NbrNodesQ1);

  /// Artificially upgrade the mesh data to higher geomtrical order
  void upgradeGeoOrder(Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes, Common::SafePtr<Common::ConnectivityTable<CFuint> > cellNodesConn, CFuint NewGeoOrder, std::vector< RealVector >& upgradedNodes,  std::vector<std::vector<cgsize_t>>& upgradedConnectivity);

  /// Clean Data by removing potential unused nodes afte Geo Upgrade
  void cleanData(std::vector<RealVector>& nodes, std::vector<std::vector<cgsize_t>>& connectivity) ;

  /// Checking for duplicate nodes
  bool nodesAreClose(const RealVector& a, const RealVector& b, int dim);
  
private:

  /// File format to write in
  std::string _fileFormatStr;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder >  m_stdTrsGeoBuilder;
  
  /// CGNS File index
  int fileIndex = -1; // Initialize to invalid value
  
  /// CGNS Base index
  int baseIndex = -1; // Initialize to invalid value  

  /// Index for CGNS zone
  int index_zone = -1;

  /// Variable to store the index of the element section in the CGNS file
  int index_section = -1;

  /// Variable to store the index of the coordinate array in the CGNS file
  int index_coord = -1;

  /// Variable to store the index of the sol array in the CGNS file
  int index_sol = -1;

  /// Dummy variable to satisfy the cg_field_write API
  int fieldIndex = -1;

  /// Helper
  std::string extractBeforeCgns(const std::string& filePath);

}; // class CGNSHighOrderWriter

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CGNSWriter_CGNSHighOrderWriter_hh

