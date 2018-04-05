// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BaseGeometricEntityProvider_hh
#define COOLFluiD_Framework_BaseGeometricEntityProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NamedObject.hh"
#include "Environment/Provider.hh"
#include "Framework/InterpolatorProperties.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
        
//////////////////////////////////////////////////////////////////////////////

/// This class represents a Provider of GeometricEntity's.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API BaseGeometricEntityProvider : public Environment::Provider<GeometricEntity>
{
 public:

  /// Default constructor without arguments
  BaseGeometricEntityProvider(const std::string& name);
  
  /// Default destructor
  virtual ~BaseGeometricEntityProvider();
  
  /// Create new geometric entity with the corresponding shape function
  /// @param states  list of the states to be put in the new GeometricEntity
  /// @param nodes   list of the nodes to be put in the new GeometricEntity
  virtual GeometricEntity* create() = 0;
  
  /// Free an instance created by this factory
  /// @param ptr pointer to be freed
  virtual void freeInstance ( void * ptr ) = 0;

  /// @return the name of this provider
  virtual std::string getProviderName () const;

  /// @return the type of this provider
  virtual std::string getProviderType () const;

  /// Get the inheritant dimensionality of the Face
  virtual CFuint getDimensionality() const = 0;

  /// Gets the type of CFGeoShape::Type
  virtual CFGeoShape::Type getShape() const = 0;

  /// Gets the type of geometry shape function
  virtual CFGeoEnt::Type getGeomType() const = 0;

  /// Gets the name of geometry shape function
  virtual std::string getGeometryShapeFunctionName() const = 0;

  /// Gets the type of geometry shape function
  virtual CFPolyForm::Type getGeometryShapeFunctionType() const = 0;

  /// Gets the geometry shape function order
  virtual CFPolyOrder::Type getGeometryShapeFunctionOrder() const = 0;

  /// Gets the name of solution shape function
  virtual std::string getSolutionShapeFunctionName() const = 0;

  /// Gets the type of solution shape function
  virtual CFPolyForm::Type getSolutionShapeFunctionType() const = 0;

  /// Gets the solution shape function  order
  virtual CFPolyOrder::Type getSolutionShapeFunctionOrder() const = 0;

  /// Gets the number of nodes associated to the geometric shape function
  virtual CFuint getGeometryShapeFunctionNbNodes() const = 0;

  /// Gets the number of nodes associated to the solution shape function
  virtual CFuint getSolutionShapeFunctionNbNodes() const = 0;

  /// Sets the ID for the solution interpolation
  virtual void setSolInterpolatorID(const InterpolatorID& id) = 0;

  /// Sets the ID for the geometric interpolation
  virtual void setGeomInterpolatorID(const InterpolatorID& id) = 0;

}; // end of class BaseGeometricEntityProvider

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BaseGeometricEntityProvider_hh
