// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IntegratorPattern_hh
#define COOLFluiD_Framework_IntegratorPattern_hh

//////////////////////////////////////////////////////////////////////////////

#include <ostream>
#include <valarray>

#include "Common/COOLFluiD.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the pattern of integration.
/// Basically it separates integration in steps,
/// each with a given number of integration points.
/// An example is a ContourIntegrator where the number of steps
/// is the number of faces in the Cell and each number per step
/// is the number of integration points in that Face.
/// @author Tiago Quintino
class Framework_API IntegratorPattern : public std::valarray<CFuint> {
public:

  /// Constructor with nb of steps as parameter
  IntegratorPattern() :  std::valarray<CFuint>(), _nbShapeF(0)
  {
  }

  /// Constructor with nb of steps as parameter
  IntegratorPattern(const CFuint& steps) :  std::valarray<CFuint>(steps), _nbShapeF(0)
  {
  }

  CFuint nbSteps() const
  {
    return size();
  }

  CFuint nbPts(const CFuint& step) const
  {
    return operator[](step);
  }

  CFuint totalNbPts() const
  {
    return sum();
  }

  std::ostream& print (std::ostream& out) const
  {
    out << "[" << *this << "]\n";
    return out;
  }

  /// Overloading of the stream operator "<<" for the output
  /// No "\n"ine introduced.
  /// @param out output stream
  /// @param p pattern to print
  /// @return the output stream
friend std::ostream& operator<< (std::ostream& out, const IntegratorPattern& p)
{
  for(CFuint i = 0; i < p.nbSteps(); ++i) {
    out << p.nbPts(i) << " ";
  }
  return out;
}

  /// @return the number of shape functions involded in the interpolation
CFuint getNbShapeFunctions() const
{
  return _nbShapeF;
}

  /// Sets the number of shape functions involded in the interpolation
  /// @param
void setNbShapeFunctions(const CFuint nbSF)
{
  _nbShapeF = nbSF;
}

private:

  /// number of Shape Functions involded in the interpolation
  CFuint _nbShapeF;

}; // end of class IntegratorPattern

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IntegratorPattern_hh
