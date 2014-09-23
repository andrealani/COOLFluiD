// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_NewtonIteratorCoupling_hh
#define COOLFluiD_Numerics_NewtonMethod_NewtonIteratorCoupling_hh

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

/// This class defines an implicit time stepper (Backward Euler like)
/// with Newton method to be used to solve weakly coupled systems of
/// equations
/// @author Andrea Lani
class NewtonIteratorCoupling : public NewtonIterator {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit NewtonIteratorCoupling(const std::string& name);

  /// Default destructor
  virtual ~NewtonIteratorCoupling();

protected:

  /// Take one timestep
  /// @see ConvergenceMethod::takeStep()
  virtual void takeStepImpl();

}; // class NewtonIteratorCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_NewtonIteratorCoupling_hh
