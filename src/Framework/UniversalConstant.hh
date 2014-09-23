// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_UniversalConstant_hh
#define COOLFluiD_Framework_UniversalConstant_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonInstantiable.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class implements the universal constants.
/// @author Steven Van Damme
/// @author Tiago Quintino
class UniversalConstant : public Common::NonInstantiable<UniversalConstant> {
public:

  static CFreal avogadro ()  { return 6.0221415e23; }
  static CFreal boltzmann () { return 1.3806503e-23; }
  static CFreal elementaryCharge () { return 1.60217646e-19; }
  static CFreal gasConstant () { return 8.31447146861745; }
  static CFreal faraday () { return 96485.3335008909; }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_UniversalConstant_hh

