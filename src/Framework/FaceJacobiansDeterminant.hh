// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FaceJacobiansDeterminant_hh
#define COOLFluiD_Framework_FaceJacobiansDeterminant_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Implements the functions for face jacobian determinant computation.
/// @author Tiago Quintino
class Framework_API FaceJacobiansDeterminant {
public:

  static CFreal compute2DFaceJacobDet(const CFreal& ds1x,
                                            const CFreal& ds1y)
  {
    return sqrt(ds1x*ds1x + ds1y*ds1y);
  }

  static CFreal compute2DFaceJacobDet(const CFreal& ds1x,
                                            const CFreal& ds1y,
                                            const CFreal& ds1z)
  {
    return sqrt(ds1x*ds1x + ds1y*ds1y + ds1z*ds1z);
  }

  static CFreal compute3DFaceJacobDet(const CFreal& ds1x,
                                            const CFreal& ds1y,
                                            const CFreal& ds1z,
                                            const CFreal& ds2x,
                                            const CFreal& ds2y,
                                            const CFreal& ds2z)
  {
    const CFreal a = (ds1y*ds2z - ds1z*ds2y);
    const CFreal b = (ds1z*ds2x - ds1x*ds2z);
    const CFreal c = (ds1x*ds2y - ds1y*ds2x);
    return sqrt(a*a + b*b + c*c);
  }

private:

  ///  Constructor
  FaceJacobiansDeterminant();

  ///  Destructor
  ~FaceJacobiansDeterminant();

}; // end class FaceJacobiansDeterminant

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FaceJacobiansDeterminant_hh
