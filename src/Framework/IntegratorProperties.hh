// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IntegratorProperties_hh
#define COOLFluiD_Framework_IntegratorProperties_hh

//////////////////////////////////////////////////////////////////////////////

#include <ostream>

#include "Common/StringOps.hh"

#include "Framework/CFPolyForm.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/CFQuadrature.hh"
#include "Framework/CFIntegration.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// IntegratorID is an ID that helps to access integrators in an unique fast way
/// without making use of virtual functions.
typedef CFuint IntegratorID;

/// This class assembles the properties of an integrator, i.e. a ShapeFunction
class Framework_API IntegratorProperties {
public:

  IntegratorProperties(const std::string&          integratName,
                       const CFIntegration::Type&  integratType,
                       const CFQuadrature::Type&   quadratureType,
                       const CFPolyOrder::Type&    integratOrder,
                       const CFGeoShape::Type&     interpolShape,
                       const CFPolyForm::Type&     interpolType,
                       const CFPolyOrder::Type&    interpolOrder);

  /// Overloading of the stream operator "<<" for the output.
  /// A "\n" is introduced at the end of every line with a property.
  friend std::ostream& operator<< (std::ostream& out, const IntegratorProperties& prop);

public:

  /// name of the integrator
  const std::string Name;
  /// integration type of the integrator
  const CFIntegration::Type IntegratType;
  /// quadrature type of the integration
  const CFQuadrature::Type QuadType;
  /// polynomial type of the interpolator
  const CFPolyForm::Type InterpolType;
  /// order of the integrator
  const CFPolyOrder::Type Order;
  /// shape on which it operates
  const CFGeoShape::Type Shape;
  /// order of the interpolator
  const CFPolyOrder::Type InterpolOrder;

}; // end class IntegratorProperties

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IntegratorID_hh
