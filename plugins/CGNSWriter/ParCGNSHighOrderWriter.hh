// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CGNSWriter_ParCGNSHighOrderWriter_hh
#define COOLFluiD_IO_CGNSWriter_ParCGNSHighOrderWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ParFileWriter.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

#include "CGNSWriter/CGWriterData.hh"

#include <unordered_set>


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

class CGNSWriter_API ParCGNSHighOrderWriter : public CGWriterCom,
                      public Framework::ParFileWriter {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit ParCGNSHighOrderWriter(const std::string& name);

  /// Destructor.
  ~ParCGNSHighOrderWriter()
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
  void openFile(const std::string& fileName, MPI_Comm comm, MPI_Info info);

  /// Write CGNS Base
  void writeCGNSBase(const CFuint cellDim, const CFuint physDim, std::string BaseName, int& fileIndex, int& baseIndex, MPI_Comm comm);

  /// Write CGNS Zone
  void writeCGNSZone(cgsize_t size[3], std::string ZoneName, int& fileIndex, int& baseIndex, int& index_zone);

  /// Write nodal coordinates
  void writeCGNSNodeCoordinates(CFuint totNbNodes, CFuint dim, std::vector< RealVector > nodes,int& fileIndex, int& baseIndex, int& index_zone, int& index_coord_x, int& index_coord_y, int& index_coord_z, cgsize_t start, cgsize_t end);

  /// Write CGNS solution node
  void writeCGNSSolutionNode(const std::string& solutionName, int& fileIndex, int& baseIndex, int& index_zone, int& index_sol);

  /// Write CGNS cell-centered solution node (for P0)
  void writeCGNSSolutionNodeCellCenter(const std::string& solutionName, int& fileIndex, int& baseIndex, int& index_zone, int& index_sol);

  /// returns the local subcell node connectivity in a cell of given shape and order (in CGNS conventions)
  static std::vector<CFuint> getCGNSCellNodeConn(ElementType_t elementType);

  /// returns the CGNS cell shape
  static ElementType_t getCGNSCellShape(CFGeoShape::Type shape,CFuint geoOrder);

  /// returns the nodal coords in the reference domain depending on geoOrder in CGNS conventions
  std::vector< RealVector > getOutputPntsMappedCoords(ElementType_t type);

  /// returns the nodal coords of the additional nodes in the physical domain
  std::vector< RealVector > calculateNewNodePositionsForElement(std::vector< RealVector >  outputPntsMappedCoords, CFuint nbrNodes, std::vector<COOLFluiD::Framework::Node*>*& cellNodes, std::vector< RealVector > geoShapeFuncs, CFuint NbrNodesQ1);

  /// Artificially upgrade the mesh data to higher geomtrical order
  void upgradeGeoOrder(std::vector<Framework::Node*> nodes, std::vector<std::vector<CFuint>> cellNodesConn, CFuint NewGeoOrder, std::vector< RealVector >& upgradedNodes,  std::vector<std::vector<cgsize_t>>& upgradedConnectivity);

  /// Clean Data by removing potential unused nodes afte Geo Upgrade
  void cleanData(std::vector<RealVector>& nodes, std::vector<std::vector<cgsize_t>>& connectivity) ;

  /// Checking for duplicate nodes
  bool nodesAreClose(const RealVector& a, const RealVector& b, int dim);

  /// Clean Nodes and Connectivity after re;oving ghost elements
  void cleanGhostData(std::vector<Framework::Node*>& nodes, std::vector<std::vector<CFuint>>& connectivity);

  /// Write P0 cell-centered solution data
  void writeP0CellCenteredSolution(MPI_Comm comm, int rank, int fileIndex, int baseIndex, 
                                    int index_zone, int index_sol, CFuint totNbCells, 
                                    CFuint globalNbCells, CFuint nbEqs,
                                    const std::vector<std::string>& varNames,
                                    const std::vector<std::string>& extraVarNames,
                                    const std::vector<std::string>& dh_varnames,
                                    Common::SafePtr<Framework::ConvectiveVarSet> updateVarSet,
                                    Common::SafePtr<Framework::DataHandleOutput> datahandle_output,
                                    Common::SafePtr<Framework::TopologicalRegionSet> trs,
                                    Framework::StdTrsGeoBuilder::GeoData& geoData,
                                    Common::SafePtr<std::vector<Framework::ElementTypeData>> elemType,
                                    CFuint nbrElemTypes);

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
  int index_coord_x = -1;
  int index_coord_y = -1;
  int index_coord_z = -1;

  /// Variable to store the index of the sol array in the CGNS file
  int index_sol = -1;

  /// Dummy variable to satisfy the cg_field_write API
  int fieldIndex = -1;

  /// Helper
  std::string extractBeforeCgns(const std::string& filePath);

}; // class ParCGNSHighOrderWriter

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CGNSWriter_ParCGNSHighOrderWriter_hh

