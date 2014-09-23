// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_InterpolatorProperties_hh
#define COOLFluiD_Framework_InterpolatorProperties_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"

#include "Framework/CFGeoEnt.hh"
#include "Framework/CFPolyForm.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// InterpolatorID is an ID that helps to access interpolators in an unique fast way
/// without making use of virtual functions.
typedef CFuint InterpolatorID;

/// This class assembles the properties of an interpolator, i.e. a ShapeFunction
class Framework_API InterpolatorProperties {
public:

  InterpolatorProperties(const std::string&    interpolName,
                         const CFPolyForm::Type&  interpolType,
                         const CFPolyOrder::Type& interpolOrder,
                         const CFGeoShape::Type&  interpolShape) :
    Name(interpolName),
    Type(interpolType),
    Order(interpolOrder),
    Shape(interpolShape)
  {
  }

public:

  /// name of the interpolator
  const std::string Name;
  /// polunomial type of the interpolator
  const CFPolyForm::Type Type;
  /// order of the interpolator
  const CFPolyOrder::Type Order;
  /// shape which it interpolates
  const CFGeoShape::Type Shape;

}; // end class InterpolatorProperties

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_InterpolatorProperties_hh
