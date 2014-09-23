// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CFmeshFileWriterFluctSplitP1P2_hh
#define COOLFluiD_Framework_CFmeshFileWriterFluctSplitP1P2_hh

//////////////////////////////////////////////////////////////////////////////



#include "FileWriter.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CFmesh format writer.
/// @author Tiago Quintino
/// @author Andrea Lani
template <class DATA>
class CFmeshFileWriterFluctSplitP1P2 : public FileWriter {
public:

  /// Constructor.
  explicit CFmeshFileWriterFluctSplitP1P2();

  /// Destructor.
  virtual ~CFmeshFileWriterFluctSplitP1P2();

  /// Sets the pointer to the stored data
  void setWriteData(const Common::SafePtr<DATA>& data)
  {
    m_writeData = data;
  }

  /// Gets the objedct where data is stored
  DATA& getWriteData()
  {
    cf_assert(m_writeData.isNotNull());
    return * m_writeData;
  }

  /// Gets the objedct where data is stored
  const DATA& getWriteData() const
  {
    cf_assert(m_writeData.isNotNull());
    return * m_writeData;
  }

  /// Get the file extension
  const std::string getWriterFileExtension() const
  {
    return std::string(".CFmesh");
  }

  /// Releases all temporary memory created while writing the file
  void releaseTemporaryWriteMemory()
  {
    getWriteData().releaseMemory();
  }

protected: // methods

  /// Writes to the given file.
  /// @throw Common::FilesystemException
  void writeToFileStream(std::ofstream& fout);

  /// Get the name of the reader
  const std::string getWriterName() const
  {
    return "CFmeshFileWriter";
  }

private: // helper functions

  /// Writes the space dimension
  void writeDimension(std::ofstream& fout);

  /// Writes the version stamp
  void writeVersionStamp(std::ofstream& fout);

  /// Writes the number of equations
  void writeNbEquations(std::ofstream& fout);

  /// Writes the extra variables info
  void writeExtraVarsInfo(std::ofstream& fout);

  /// Writes the number of nodes
  void writeNbNodes(std::ofstream& fout);

  /// Writes the nb of dofs state tensors and initialize with them the dofs
  void writeNbStates(std::ofstream& fout);

  /// Writes the nb of elements
  void writeNbElements(std::ofstream& fout);

  /// Writes the nb of element types
  void writeNbElementTypes(std::ofstream& fout);

  /// Writes the order of the polynomial representation of the geometry
  void writeGeometricPolyOrder(std::ofstream& fout);

  /// Writes the order of the polynomial representation of the solution
  void writeSolutionPolyOrder(std::ofstream& fout);

  /// Writes the element types (CFGeoShape::Type)
  void writeElementTypes(std::ofstream& fout);

  /// Writes the nb of elements per type
  void writeNbElementsPerType(std::ofstream& fout);

  /// Writes the nb of nodes per type
  void writeNbNodesPerType(std::ofstream& fout);

  /// Writes the nb of dofs per type
  void writeNbStatesPerType(std::ofstream& fout);

  /// Writes the list of nodes
  void writeNodeList(std::ofstream& fout);

  /// Writes the list of state tensors and initialize the dofs
  void writeStateList(std::ofstream& fout);

  /// Writes the data concerning the elements
  void writeElementList(std::ofstream& fout);

  /// Writes the number of topological region sets and initialize the vector
  /// that will contain the all the topological region sets
  void writeNbTRSs(std::ofstream& fout);

  /// Writes the all the data relative to all TRSs
  void writeTrsData(std::ofstream& fout);

protected: // data

  /// acquaintance of the data present in the CFmesh file
  Common::SafePtr<DATA> m_writeData;

}; // class CFmeshFileWriterFluctSplitP1P2

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileWriterFluctSplitP1P2.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CFmeshFileWriterFluctSplitP1P2_hh

