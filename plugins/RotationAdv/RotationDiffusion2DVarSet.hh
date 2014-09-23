// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdv_RotationDiffusion2DVarSet_hh
#define COOLFluiD_Physics_RotationAdv_RotationDiffusion2DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "RotationDiffusionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables
   *
   * @author Nadege Villedieu
   */
class RotationDiffusion2DVarSet : public RotationDiffusionVarSet {
public: // classes

  /**
   * Constructor
   * @see RotationDiffusion2D
   */
  RotationDiffusion2DVarSet(const std::string& name,
		       Common::SafePtr<Framework::PhysicalModelImpl> model) :
    RotationDiffusionVarSet(name, model)
  {
  }

  /**
   * Default destructor
   */
  virtual ~RotationDiffusion2DVarSet()
  {
  }

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);
  
  
}; // end of class RotationDiffusion2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_RotationAdv_RotationDiffusion2DVarSet_hh
