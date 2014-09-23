// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_NewmarkExplicit_hh
#define COOLFluiD_Numerics_NewtonMethod_NewmarkExplicit_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the NewmarkExplicit
/// method for solving non linear systems
/// @author Thomas Wuilbaut
class NewmarkExplicit : public NewtonIterator {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit NewmarkExplicit(const std::string& name);

  /// Default destructor
  ~NewmarkExplicit();

protected:

  /// Take one timestep
  /// @see ConvergenceMethod::takeStep()
  virtual void takeStepImpl();

}; // class NewmarkExplicit

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewmarkExplicit_hh
