// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_LinearAdvSys3DLinearPrim_hh
#define COOLFluiD_Physics_LinearAdvSys_LinearAdvSys3DLinearPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/JacobianLinearizer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class State; }

  namespace Physics {

    namespace LinearAdvSys {

      class LinearAdvSysTerm;

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class performs all the operations related with physics-dependent
  * transformations of variables
  *
  * @author Nadege Villedieu
  */
class LinearAdvSys3DLinearPrim : public Framework::JacobianLinearizer {
public:

  /**
   * Default constructor without arguments
   */
  LinearAdvSys3DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model);

  /**
   * Default destructor
   */
  ~LinearAdvSys3DLinearPrim();

  /**
   * Linearize the states
   */
  void linearize(const std::vector<Framework::State*>& statesInCell);

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<LinearAdvSysTerm> _model;

}; // end of class LinearAdvSys3DLinearPrim

//////////////////////////////////////////////////////////////////////////////

    } //  namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Euler_Euler3DLinearPrim_hh
