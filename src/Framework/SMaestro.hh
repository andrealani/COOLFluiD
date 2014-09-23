// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SMaestro_hh
#define COOLFluiD_Framework_SMaestro_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Maestro.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a sequential Maestro.
/// It sequentially solves each SubSystem once (or more).
/// It controls the flow of actions in the SubSystem's through Event's.
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API SMaestro : public Framework::Maestro {
public:

  /// Default constructor without arguments.
  SMaestro(const std::string& name);

  /// Default destructor.
  ~SMaestro();

  /// Takes control of the simulation
  Common::Signal::return_t control ( Common::Signal::arg_t );

}; // end of class SMaestro

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SMaestro_hh
