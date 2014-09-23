// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CFmeshFileReader_hh
#define COOLFluiD_Framework_CFmeshFileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "Common/SafePtr.hh"
#include "Common/Table.hh"

#include "MathTools/RealVector.hh"

#include "Framework/FileReader.hh"
#include "Framework/ElementTypeData.hh"
#include "Framework/BadFormatException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CFmesh format reader.
/// @author Tiago Quintino
/// @author Andrea Lani
template <class DATA>
class CFmeshFileReader : public FileReader {

public: // functions

  /// Constructor.
  CFmeshFileReader();

  /// Destructor.
  ~CFmeshFileReader();

  /// Sets the pointer to the stored data
  void setReadData(const Common::SafePtr<DATA>& data)
  {
    m_readData = data;
  }

  /// Gets the object where data is stored
  DATA& getReadData()
  {
    cf_assert(m_readData.isNotNull());
    return *m_readData;
  }

  /// Gets the objedct where data is stored
  const DATA& getReadData() const
  {
    cf_assert(m_readData.isNotNull());
    return *m_readData;
  }

  /// Get the file extension
  const std::string getReaderFileExtension() const
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

private: // typedefs

  /// pointer to ReaderFunction
  typedef void (CFmeshFileReader<DATA>::*ReaderFunction)(std::ifstream& fin);

  /// type that maps a string read in the File with a ReaderFunction
  typedef std::map<std::string,
      ReaderFunction,
      std::less<std::string> > MapString2Reader;

private: // helper functions

  /// Read an entry in the .CFmesh file
  bool readString(std::ifstream& file);

  /// Get the file extension
  const std::string getReaderTerminator() const
  {
    static const std::string terminator = "!END";
    return terminator;
  }

  /// Get the name of the reader
  const std::string getReaderName() const
  {
    return "CFmeshFileReader";
  }

  /// Sets m_mapString2Reader, that maps a given string to a corresponding
  /// reader function
  void setMapString2Readers();

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

  /// Reads the nb of extra variables associated with nodes
  void readNbExtraNodalVars(std::ifstream& fin);

  /// Reads the nb of extra variables associated with states
  void readNbExtraStateVars(std::ifstream& fin);

  /// Reads the names of the extra variables associated with nodes
  void readExtraNodalVarNames(std::ifstream& fin);

  /// Reads the names of the extra variables associated with states
  void readExtraStateVarNames(std::ifstream& fin);

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

private: // data

  /// acquaintance of the data present in the CFmesh file
  Common::SafePtr<DATA>       m_readData;

  /// map each string with a corresponding pointer to member
  /// function
  MapString2Reader           m_mapString2Reader;

  /// array of flags telling if to use the given initial
  /// values to initialize each variable in the State's
  std::vector<bool> m_useInitValues;

  /// array of initial values to set in the State's
  std::vector<CFreal> m_initValues;

  /// array of initial values IDs to set in the State's
  std::vector<CFuint> m_initValuesIDs;

}; // class CFmeshFileReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileReader.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CFmeshFileReader_hh

