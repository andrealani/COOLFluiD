// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_LinearAdv3DLinearPrim_hh
#define COOLFluiD_Physics_LinearAdv_LinearAdv3DLinearPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {   class State;  }

  namespace Physics {
  namespace LinearAdv {

    class LinearAdvTerm;

//////////////////////////////////////////////////////////////////////////////

/// This class calculates the useful physical-model dependent
/// quantities in a linearized state
/// @author Tiago Quintino
class LinearAdv_API LinearAdv3DLinearPrim : public Framework::JacobianLinearizer {
public:

  /// Default constructor without arguments
  LinearAdv3DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model);

  /// Default destructor
  ~LinearAdv3DLinearPrim();

  /// Linearize the states
  void linearize(const std::vector<Framework::State*>& statesInCell);

private:

  /// physical model
  Common::SafePtr<LinearAdvTerm> _model;

}; // end of class LinearAdv3DLinearPrim

//////////////////////////////////////////////////////////////////////////////

    } //  namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_LinearAdv3DLinearPrim_hh
