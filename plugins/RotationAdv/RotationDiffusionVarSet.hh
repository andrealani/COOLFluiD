// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdv_RotationDiffusionVarSet_hh
#define COOLFluiD_Physics_RotationAdv_RotationDiffusionVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "RDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model  for primitive
   * variables
   *
   * @author Nadege Villedieu
   */
class RotationDiffusionVarSet : public Framework::DiffusiveVarSet {
public: // classes

  /**
   * Constructor
   * @see RotationDiffusion
   */
  RotationDiffusionVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model->getDiffusiveTerm().d_castTo<RDTerm>()),
    m_diffcoeff(0.0)
  {
  }

  /**
   * Default destructor
   */
  virtual ~RotationDiffusionVarSet()
  {
  }

  /**
   * Set up private data
   */
  virtual void setup()
  {
    Framework::DiffusiveVarSet::setup();
    m_gradState.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  }

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  virtual void setComposition(const RealVector& state,
                              const bool isPerturb,
                              const CFuint iVar) = 0;

  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  virtual void setGradientVars(const std::vector<RealVector*>& states,
                               const std::vector<RealVector*>& values,
                               const CFuint stateSize) = 0;

    /**
   * Get the current diffusion coefficient
   */
  virtual CFreal getCurrDiffuionCoeff()
  {
    return (getModel().getPhysicalData())[RDTerm::NU];
  }

  /**
   * Get the model
   */
  RDTerm& getModel()
  {
    return *_model;
  }

  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius) = 0;

  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius) = 0;

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state) = 0;

private:

  /// physical model
  Common::SafePtr<RDTerm> _model;

protected:

  /// gradient variables
  RealVector m_gradState;

  /// Diffusion coefficient
  CFreal m_diffcoeff;

}; // end of class RotationDiffusionVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_RotationAdv_RotationDiffusionVarSet_hh
