// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshDataBuilder_hh
#define COOLFluiD_Framework_MeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Common/CFMultiMap.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"

#include "Config/ConfigObject.hh"

#include "Environment/ConcreteProvider.hh"

#include "Framework/CFPolyForm.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFmeshReaderSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

  class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

/// This class offers an abstract interface for MeshDataBuilder's.
/// ReadCFmesh aggregates the MeshDataBuilder.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API MeshDataBuilder : public Common::OwnedObject,
				      public Common::NonCopyable<MeshDataBuilder>,
				      public Config::ConfigObject
{

public: // typedefs

  typedef Environment::ConcreteProvider<MeshDataBuilder,1> PROVIDER;
  typedef const std::string& ARG1;

public: // static functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

public: // functions

  /// Constructor
  /// @param name name of the builder used for configuration
  MeshDataBuilder(const std::string& name);

  /// Virtual destructor
  virtual ~MeshDataBuilder();

  /// Configures this object from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Creates all the TopologicalRegionSet's from CFmeshData
  virtual void createTopologicalRegionSets();

  /// Create the inner cells.
  /// This function provides a default implementation.
  virtual void createInnerCells();

  /// Creates the TRS of the boundaries.
  /// This function provides a default implementation.
  virtual void createBoundaryTRS();

  /// Creates a the TopologicalRegionSet and put it in the MeshDataBuilder
  /// @param nbGeomEntsInTr  list of the nb of geometric entities in each tr
  /// @param name  name of the TRS
  /// @param trGeoCon connectivity for GeometricEntity's present in the TRS
  virtual Common::SafePtr<Framework::TopologicalRegionSet> createTopologicalRegionSet
  (const std::vector<CFuint>& nbGeomEntsInTr,
  const std::string& name,
  const TRGeoConn& trGeoCon,
  const CFuint iTRS);

  /// Compute the all the needed connectivity data
  /// @pre setData() has to be called before this method
  /// Compute the connectivity
  virtual void computeGeoTypeInfo();

  /// Releases temporary memory used in building the mesh
  /// Should be called by child classes
  virtual void releaseMemory();

  /// Set some global info like the maximum
  /// number of dofs or faces in the cells
  void setMaxGlobalInfo()
  {
    setMaxNbStatesInCell();
    setMaxNbNodesInCell();
    setMaxNbFacesInCell();
  }

  /// Set CFmeshData needed by the MeshDataBuilder
  void setCFmeshData(Common::SafePtr<CFmeshReaderSource> data);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "MeshDataBuilder";
  }

protected: // methods

  /// Set the max number of states in cell
  virtual void setMaxNbStatesInCell() = 0;

  /// Set the max number of nodes in cell
  virtual void setMaxNbNodesInCell() = 0;

    /// Set the max number of faces in cell
  virtual void setMaxNbFacesInCell() = 0;

protected: // helper methods

  /// @return the number of elements
  CFuint getNbElements() const;

  /// Get the number of boundary faces
  CFuint getNbBoundaryFaces() const;

  /// @return the number of element types
  CFuint getNbElementTypes() const;

  /// @return the number of nodes
  CFuint getNbNodes() const;

  /// @return the number of states
  CFuint getNbStates() const;

  /// @return the order of the geometrical polynomial representation
  /// of the mesh
  CFPolyOrder::Type getGeometricPolyOrder() const;

  /// @return the order of the polynomial representation of the solution
  CFPolyOrder::Type getSolutionPolyOrder() const;

  /// @return the polynomial type for the geometry
  CFPolyForm::Type getGeometricPolyType() const;

  /// @return the polynomial type for the solution
  CFPolyForm::Type getSolutionPolyType() const;

  /// @return the CFGeoShape::Type corresponding to the given number of nodes
  CFGeoShape::Type getGeoShape (CFGeoEnt::Type  geomType,
    const CFuint&       dim,
                          const CFuint& nbGeoNodes) const;

  /// @return a string corresponding to the name of the GeometricEntityProvider
  std::string makeGeomEntName (const CFGeoShape::Type& shape,
                            const CFPolyForm::Type&  geomPolyType,
                            const CFPolyOrder::Type& geomPolyOrder,
                            const CFPolyForm::Type&  solPolyType,
                            const CFPolyOrder::Type& solPolyOrder) const;

  ///  Gets the CFmeshData used in the building of MeshData
  ///  (write access)
  CFmeshReaderSource& getCFmeshData()
  {
    cf_assert(m_cfmeshData.isNotNull());
    return *m_cfmeshData;
  }

  ///  Gets the CFmeshData used in the building of MeshData
  /// (ONLY read access)
  const CFmeshReaderSource& getCFmeshData() const
  {
    cf_assert(m_cfmeshData.isNotNull());
    return *m_cfmeshData;
  }

  ///  Set the element type data in the MeshData
  void setElementTypeData();

  ///  Set the coordinates (Node's) in cell State's
  void setCoordInCellStates();

  /// Set the arrays storing the number of states and nodes for
  /// all GeometricEntity's in the current TRS
  void setNbNodesNbStatesInGeo(const std::vector<CFuint>& nbGeomEntsPerTR,
			       const TRGeoConn& trGeoConn,
			       std::valarray<CFuint>& nbNodesInGeo,
			       std::valarray<CFuint>& nbStatesInGeo);
  
private: // methods

  /// Create and set the mapping between faces and TRSs.
  /// Allows to get Face info from TRS by geometric entity ID,
  /// local to the processor.
  virtual void setMapGeoToTrs() = 0;

private: // data

  /// acquaintance with the CFmesh data
  Common::SafePtr<CFmeshReaderSource> m_cfmeshData;

  /// polynomial type of the geometric entities to be built
  std::string m_polynomialTypeName;

protected: // data

  /// current boundary face ID
  CFuint m_currBFaceID;

  /// map between the GeometricEntity's provider name and the
  /// corresponding type
  Common::CFMap<std::string, CFuint> m_mapGeoProviderNameToType;

}; // end of class MeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(MeshDataBuilder) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshDataBuilder_hh
