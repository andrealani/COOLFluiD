// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_AdvectionDiffusionSys3DPrim_hh
#define COOLFluiD_Physics_LinearAdvSys_AdvectionDiffusionSys3DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "AdvectionDiffusionSys3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdvSys {

      class LinearAdvSysTerm;

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NavierStokes physical model 2D for primitive
  /// variables
  /// @author Nadege Villedieu
class AdvectionDiffusionSys3DPrim : public AdvectionDiffusionSys3DVarSet {
public: // classes

  /// Constructor
  /// @see AdvectionDiffusionSys2D
  AdvectionDiffusionSys3DPrim(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  ~AdvectionDiffusionSys3DPrim();

  /// Set the composition
  /// @pre this function has to be called before any other function
  ///      computing other physical quantities
  void setComposition(const RealVector& state,
  	      const bool isPerturb,
  	      const CFuint iVar)
  {
  }

  /// Set the quantities needed to compute gradients (pressure,
  /// velocity, etc.) starting from the states
  void setGradientVars(const std::vector<RealVector*>& states,
					 RealMatrix& values,
		       const CFuint stateSize);


protected:

  /// Set the gradient variables starting from state variables
  virtual void setGradientState(const RealVector& state);

private:

  /// convective model
  Common::SafePtr<LinearAdvSysTerm> _linearadvModel;

}; // end of class AdvectionDiffusion3DPrim

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_AdvectionDiffusion3DPrim_hh
