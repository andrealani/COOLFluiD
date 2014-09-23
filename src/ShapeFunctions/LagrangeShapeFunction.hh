// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_LagrangeShapeFunction_hh
#define COOLFluiD_ShapeFunctions_LagrangeShapeFunction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ShouldNotBeHereException.hh"
#include "Common/NotImplementedException.hh"

#include "MathTools/RealVector.hh"

#include "Framework/Node.hh"
#include "Framework/FaceJacobiansDeterminant.hh"
#include "Framework/IntegratorPattern.hh"
#include "Framework/InterpolatorProperties.hh"
#include "Framework/PhysicalModel.hh"

#include "ShapeFunctions/ShapeFunctionsAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

/// This class provides methods and data common to all the
/// lagrangian shape functions
/// @author Andrea Lani
/// @author Tiago Quintino
class ShapeFunctions_API LagrangeShapeFunction {
public:

  /// Interpolate the solution.
  /// This is a Template Method.
  /// @param nodes   list of Node's or State's in a given GeometricEntity
  /// @param nodalSF nodal shape functions in the given GeometricEntity
  /// @param values  array in which to put the interpolated values
  template <class T>
  static inline void interpolate(const std::vector<T*>& nodes,
                                 const std::vector<RealVector>& nodalSF,
                                       std::vector<T*>& values);

protected:

  /// Default constructor without arguments
  LagrangeShapeFunction();

  /// Default destructor
  ~LagrangeShapeFunction();

}; // end of class LagrangeShapeFunction

//////////////////////////////////////////////////////////////////////////////

template <class T>
inline void LagrangeShapeFunction::interpolate(
       const std::vector<T*>& nodes,
       const std::vector<RealVector>& nodalSF,
             std::vector<T*>& values)
{
  using namespace std;

  const CFuint nbInterpolationPoints = nodalSF.size();
  const CFuint nbNodes = nodes.size();
  cf_assert(nodes.size() > 0);
  cf_assert(values.size() >= nbInterpolationPoints);

  for (CFuint ip = 0; ip < nbInterpolationPoints; ++ip)
  {
    *values[ip] = 0.;
    cf_assert(nodalSF[ip].size() == nodes.size());
    for (CFuint n = 0; n < nbNodes; ++n)
    {
      cf_assert(values[ip]->size()  == nodes[n]->size());
      for (CFuint iVar = 0; iVar < nodes[n]->size(); ++iVar)
      {
        (*values[ip])[iVar] += nodalSF[ip][n]*(*nodes[n])[iVar];
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////

} // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_LagrangeShapeFunction_hh
