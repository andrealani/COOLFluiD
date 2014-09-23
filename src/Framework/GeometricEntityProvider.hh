// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GeometricEntityProvider_hh
#define COOLFluiD_Framework_GeometricEntityProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/BaseGeometricEntityProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class GeometricEntity;
    class State;
    class Node;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a factory of GeometricEntity, Cell's
/// and Face's, with ShapeFunction's
/// The class is parameterized with the type of geometric entity, e.g. Cell,
/// Face and  Edge
/// @author Andrea Lani
template <template
         <class GEO_SHAPE_FUNCTION,
          class SOL_SHAPE_FUNCTION> class GEOENTTYPE,
          class GEO_SHAPE_FUNCTION,
          class SOL_SHAPE_FUNCTION,
          class LIBRARY >
class GeometricEntityProvider : public BaseGeometricEntityProvider {
public:

  /// Constructor
  /// @param name
  /// @see GeometricEntityProvider()
  GeometricEntityProvider(const std::string& name) : BaseGeometricEntityProvider(name) 
  {
    Environment::ModuleRegister<LIBRARY>::getInstance().getSelfRegistry().regist(this);
  }

  /// Default destructor
  ~GeometricEntityProvider() {}

  /// Create new geometric entity with the corresponding shape function
  /// @param states list of the states to be put in the new GeometricEntity
  /// @param nodes  list of the nodes to be put in the new GeometricEntity
  GeometricEntity* create ()
  {
    return new GEOENTTYPE<GEO_SHAPE_FUNCTION,SOL_SHAPE_FUNCTION>();
  }

  /// Free an instance created by this provider
  /// @param ptr pointer to be freed
  void freeInstance ( void* ptr )
  {
    cf_assert(ptr != CFNULL);
    GeometricEntity* obj = reinterpret_cast<GeometricEntity*>(ptr);

    cf_assert(obj != CFNULL);
    deletePtr<GeometricEntity>( obj );
  }

  /// Get the inheritant dimensionality of the Face
  CFuint getDimensionality() const
  {
    cf_assert(GEO_SHAPE_FUNCTION::getDimensionality() == SOL_SHAPE_FUNCTION::getDimensionality());
    return GEO_SHAPE_FUNCTION::getDimensionality();
  }

  /// Gets the CFGeoEnt::Type of the provided GeometricEntity
  /// @return CFGeoEnt::Type
  CFGeoEnt::Type getGeomType() const
  {
    cf_assert(GEO_SHAPE_FUNCTION::getShape() == SOL_SHAPE_FUNCTION::getShape());
    return GEOENTTYPE<GEO_SHAPE_FUNCTION,SOL_SHAPE_FUNCTION>::getGeomType();
  }

  /// Gets the type of CFGeoShape::Type
  CFGeoShape::Type getShape() const
  {
    cf_assert(GEO_SHAPE_FUNCTION::getShape() == SOL_SHAPE_FUNCTION::getShape());
    return GEO_SHAPE_FUNCTION::getShape();
  }

  /// Gets the name of geometry shape function
  std::string getGeometryShapeFunctionName() const
  {
    return GEO_SHAPE_FUNCTION::getName();
  }

  /// Gets the type of geometry shape function
  CFPolyForm::Type getGeometryShapeFunctionType() const
  {
    return GEO_SHAPE_FUNCTION::getInterpolatorType();
  }

  /// Gets the geometry shape function order
  CFPolyOrder::Type getGeometryShapeFunctionOrder() const
  {
    return GEO_SHAPE_FUNCTION::getInterpolatorOrder();
  }

  /// Gets the name of solution shape function
  std::string getSolutionShapeFunctionName() const
  {
    return SOL_SHAPE_FUNCTION::getName();
  }

  /// Gets the type of solution shape function
  CFPolyForm::Type getSolutionShapeFunctionType() const
  {
    return SOL_SHAPE_FUNCTION::getInterpolatorType();
  }

  /// Gets the solution shape function  order
  CFPolyOrder::Type getSolutionShapeFunctionOrder() const
  {
    return SOL_SHAPE_FUNCTION::getInterpolatorOrder();
  }

  /// Gets the number of nodes associated to the geometric shape function
  CFuint getGeometryShapeFunctionNbNodes() const
  {
    return GEO_SHAPE_FUNCTION::getNbNodes();
  }

  /// Gets the number of nodes associated to the solution shape function
  CFuint getSolutionShapeFunctionNbNodes() const
  {
    return SOL_SHAPE_FUNCTION::getNbNodes();
  }

  /// Sets the ID for the solution interpolation
  void setSolInterpolatorID(const InterpolatorID& id)
  {
    SOL_SHAPE_FUNCTION::setInterpolatorID(id);
  }

  /// Sets the ID for the geometric interpolation
  void setGeomInterpolatorID(const InterpolatorID& id)
  {
    GEO_SHAPE_FUNCTION::setInterpolatorID(id);
  }

}; // end of class GeometricEntityProvider

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GeometricEntityProvider_hh
