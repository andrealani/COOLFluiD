// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_AdvectionDiffusion3DVarSet_hh
#define COOLFluiD_Physics_LinearAdv_AdvectionDiffusion3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "AdvectionDiffusionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NavierStokes physical model 2D for primitive
  /// variables
  /// @author Nadege Villedieu
class LinearAdv_API AdvectionDiffusion3DVarSet : public AdvectionDiffusionVarSet {
public: // classes

  /// Constructor
  /// @see AdvectionDiffusion3D
  AdvectionDiffusion3DVarSet(const std::string& name,
  	       Common::SafePtr<Framework::PhysicalModelImpl> model) :
    AdvectionDiffusionVarSet(name, model)
  {
  }

  /// Default destructor
  virtual ~AdvectionDiffusion3DVarSet()
  {
  }

  /// Set up private data
  virtual void setup();

  /// Get the diffusive flux
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

  /// Get the diffusive flux vector
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);


}; // end of class AdvectionDiffusion3DVarSet

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_AdvectionDiffusion3DVarSet_hh
