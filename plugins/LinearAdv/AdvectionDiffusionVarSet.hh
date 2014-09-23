// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_AdvectionDiffusionVarSet_hh
#define COOLFluiD_Physics_LinearAdv_AdvectionDiffusionVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "ADTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NavierStokes physical model  for primitive
  /// variables
  /// @author Nadege Villedieu
class LinearAdv_API AdvectionDiffusionVarSet : public Framework::DiffusiveVarSet {
public: // classes

  /// Constructor
  /// @see AdvectionDiffusion
  AdvectionDiffusionVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model->getDiffusiveTerm().d_castTo<ADTerm>()),
    m_diffcoeff(0.0)
  {
  }

  /// Default destructor
  virtual ~AdvectionDiffusionVarSet()
  {
  }

  /// Set up private data
  virtual void setup()
  {
    Framework::DiffusiveVarSet::setup();
    m_gradState.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  }

  /// Set the composition
  /// @pre this function has to be called before any other function
  ///      computing other physical quantities
  virtual void setComposition(const RealVector& state,
                              const bool isPerturb,
                              const CFuint iVar) = 0;

  /// Set the quantities needed to compute gradients (pressure,
  /// velocity, etc.) starting from the states
  virtual void setGradientVars(const std::vector<RealVector*>& states,
                               const std::vector<RealVector*>& values,
                               const CFuint stateSize) = 0;

  /// Get the current diffusion coefficient
  virtual CFreal getCurrDiffuionCoeff()
  {
    return (getModel().getPhysicalData())[ADTerm::NU];
  }

  /// Get the model
  ADTerm& getModel()
  {
    return *_model;
  }

  /// Get the diffusive flux
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius) = 0;

  /// Get the diffusive flux vector
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius) = 0;

protected:

  /// Set the gradient variables starting from state variables
  virtual void setGradientState(const RealVector& state) = 0;

private:

  /// physical model
  Common::SafePtr<ADTerm> _model;

protected:

  /// gradient variables
  RealVector m_gradState;

  /// Diffusion coefficient
  CFreal m_diffcoeff;

}; // end of class AdvectionDiffusionVarSet

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_AdvectionDiffusionVarSet_hh
