// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdv_RotationAdv2DLinearPrim_hh
#define COOLFluiD_Physics_RotationAdv_RotationAdv2DLinearPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace RotationAdv {
      class RotationAdvTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class calculates the useful physical-model dependent
 * quantities in a linearized state
 *
 * @author Tiago Quintino
 *
 */
class RotationAdv2DLinearPrim : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  RotationAdv2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~RotationAdv2DLinearPrim();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

 private:

  /// physical model
  Common::SafePtr<RotationAdvTerm> _model;

}; // end of class RotationAdv2DLinearPrim

//////////////////////////////////////////////////////////////////////////////

    } //  namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_RotationAdv_RotationAdv2DLinearPrim_hh
