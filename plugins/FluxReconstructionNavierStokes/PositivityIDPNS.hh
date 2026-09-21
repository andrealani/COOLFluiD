// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_PositivityIDPNS_hh
#define COOLFluiD_FluxReconstructionMethod_PositivityIDPNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/BasePositivityIDP.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Invariant domain preserving positivity limiter for Euler and Navier-Stokes.
 *
 * Supplies the physics hooks of BasePositivityIDP for the single species
 * layout (rho, rhoU, rhoV, [rhoW], rhoE), in 2D and 3D. Density and pressure
 * are constrained; there is nothing to freeze, so the Hydro and Full scaling
 * modes both scale every component and the Auto mode reduces to either.
 *
 * Update variables Cons and Prim are supported. Temperature based sets such
 * as Puvt, and any model carrying extra transported scalars such as
 * turbulence or species, are refused rather than silently constrained on the
 * wrong slots.
 *
 * @author Rayan Dhib
 */
class PositivityIDPNS : public BasePositivityIDP {

public: // functions

  /// Constructor
  explicit PositivityIDPNS(const std::string& name);

  /// Destructor
  virtual ~PositivityIDPNS();

  /// Defines the Config Option's of this class
  static void defineConfigOptions(Config::OptionList& options);

protected: // physics hooks

  /// density and pressure of a conservative state, B2 is always zero here
  virtual void constraintsAtPoint(const RealVector& cons,
                                  CFreal& rho, CFreal& p, CFreal& B2) const;

  /// map a conservative state back to the configured update variables
  virtual void consToUpdate(const RealVector& cons, RealVector& update) const;

  /// conservative component indices scaled in the given mode
  virtual const std::vector< CFuint >& scaledIndices(ScaleMode mode) const;

  /// read gamma, validate the variable set and the equation layout
  virtual void setupPhysics();

protected: // data

  /// gamma - 1 of the Euler term
  CFreal m_gammaMinusOne;

  /// component indices scaled in either mode, all of them for this physics
  std::vector< CFuint > m_allIndices;

  /// index of the total energy slot
  CFuint m_iE;

  /// true when the update variables are the conservative ones
  bool m_consUpdateVars;

}; // class PositivityIDPNS

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PositivityIDPNS_hh
