// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IntegratorTypes_hh
#define COOLFluiD_Framework_IntegratorTypes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class groups some data needed to create the geometric entity
/// and configure the right integrators
/// @author Andrea Lani
class Framework_API IntegratorTypes {
public:

  // tell if to set the contour integrator for solution shape function
  bool setContourIntegratorSS;

  // tell if to set the contour integrator for geometric shape function
  bool setContourIntegratorSG;

  // tell if to set the volume integrator for solution shape function
  bool setVolumeIntegratorSS;

  // tell if to set the volume integrator for geometric shape function
  bool setVolumeIntegratorSG;

}; // end class IntegratorTypes

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IntegratorTypes_hh
