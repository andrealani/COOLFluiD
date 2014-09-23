// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullLinearizer_hh
#define COOLFluiD_Framework_NullLinearizer_hh

//////////////////////////////////////////////////////////////////////////////

#include "JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
    class PhysicalModelImpl;

//////////////////////////////////////////////////////////////////////////////

/// This class performs all the operations related with physics-dependent
/// transformations of variables
/// @author Tiago Quintino

class Framework_API NullLinearizer : public JacobianLinearizer {

public:

  /// Default constructor without arguments
  NullLinearizer(Common::SafePtr<Framework::PhysicalModel> model);

  /// Default destructor
  ~NullLinearizer();

  /// Set the physical model
  void setPhysicalModel(PhysicalModelImpl* model);

  /// Linearize the states
  void linearize(const std::vector<State*>& statesInCell);

  /// Checks if this object is a Null object.
  /// Since this is a NullSplitter mit returns true.
  bool isNull() const
  {
    return true;
  }

}; // end of class NullLinearizer

//////////////////////////////////////////////////////////////////////////////

    } //  namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullLinearizer_hh
