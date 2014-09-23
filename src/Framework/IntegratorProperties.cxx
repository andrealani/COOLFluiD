// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/IntegratorProperties.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

IntegratorProperties::
IntegratorProperties(const std::string&          integratName,
                     const CFIntegration::Type& integratType,
                     const CFQuadrature::Type&  quadratureType,
                     const CFPolyOrder::Type&       integratOrder,
                     const CFGeoShape::Type&        interpolShape,
                     const CFPolyForm::Type&        interpolType,
                     const CFPolyOrder::Type&       interpolOrder) :
    Name(integratName),
    IntegratType(integratType),
    QuadType(quadratureType),
    InterpolType(interpolType),
    Order(integratOrder),
    Shape(interpolShape),
    InterpolOrder(interpolOrder)
{
}

//////////////////////////////////////////////////////////////////////////////

std::ostream&
operator<< (std::ostream& out, const IntegratorProperties& prop)
{

//    out << "Integrator name     : " << prop.Name << "\n";
//    out << "Integration type    : " << CFIntegration::Convert::to_str(prop.IntegratType) << "\n";
//    out << "Geometric shape     : " << CFGeoShape::Convert::to_str(prop.Shape) << "\n";
//    out << "Interpolation type  : " << CFPolyForm::Convert::to_str(prop.InterpolType) << "\n";
//    out << "Interpolation order : " << CFPolyOrder::Convert::to_str(prop.InterpolOrder) << "\n";
//    out << "Quadrature type     : " << CFQuadrature::Convert::to_str(prop.QuadType) << "\n";
//    out << "Integration order   : " << CFPolyOrder::Convert::to_str(prop.Order) << "\n";

   out << "[" << prop.Name  << "] ";
   out << "[" << prop.Shape << "] ";
   out << "[" << prop.IntegratType << "] ";
   out << "[" << prop.Order << "] ";
   out << "[" << prop.QuadType << "] ";
   out << "[" << prop.InterpolType << "] ";
   out << "[" << prop.InterpolOrder << "] ";

   return out;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

