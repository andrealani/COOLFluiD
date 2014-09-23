// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdv_RotationDiffusion2DPrim_hh
#define COOLFluiD_Physics_RotationAdv_RotationDiffusion2DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "RotationDiffusion2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {
      class RotationAdvTerm;

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables
   *
   * @author Nadege Villedieu
   */
class RotationDiffusion2DPrim : public RotationDiffusion2DVarSet {
public: // classes

  /**
   * Constructor
   * @see RotationDiffusion2D
   */
  RotationDiffusion2DPrim(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~RotationDiffusion2DPrim();

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  void setComposition(const RealVector& state,
		      const bool isPerturb,
		      const CFuint iVar)
  {
  }

  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  void setGradientVars(const std::vector<RealVector*>& states,
		       const std::vector<RealVector*>& values,
		       const CFuint stateSize);


protected:
  
  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

private:

  /// convective model
  Common::SafePtr<RotationAdvTerm> _RotationAdvModel;

}; // end of class RotationDiffusion2DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_RotationAdv_RotationDiffusion2DPrim_hh
