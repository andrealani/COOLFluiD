// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSolCorona_hh
#define COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSolCorona_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonMethod/PhysicalityCFLUpdateSol.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// PhysicalityCFLUpdateSol with the solar corona state clipping of
/// StdUpdateSolPP2 applied after an accepted update. The step control, the
/// CFL schedule and the retries are the ones of the base class, this command
/// only adds, per parallel updatable state and in this order:
///  - an Alfven speed cap, which raises the density wherever the local Alfven
///    speed runs past AlfvenFactor times the reference speed built from the
///    largest field strength in the domain and DensityBC,
///  - a temperature floor at fixed density, which raises the pressure,
///  - a velocity damping below DampingRadius, which allows no more than the
///    local sound speed at the inner boundary and releases it with height.
/// Each cap is approached through a tanh ramp, so a state sitting under one
/// is left alone.
///
/// The clipping only ever raises rho and p and shrinks the velocity, so it
/// cannot undo the physicality the base class established for the step.
///
/// The three treatments index the state directly as
/// [rho, u, v, w, Bx, By, Bz, p, ...], which MHD3DProjectionPrim and
/// MHD3DProjectionPrimE provide, and read the corona normalization from
/// RefValues. The variable names are checked on the first update. When a
/// potential field is active, the B entering the Alfven speed is the
/// deviation field the state carries, as in StdUpdateSolPP2.
///
/// One difference with StdUpdateSolPP2: the largest field strength that sets
/// the Alfven speed cap is reduced over all ranks here. StdUpdateSolPP2 keeps
/// it rank local, which makes its density cap depend on the partitioning.
///
/// The treatment itself is the COCONUT one, from StdUpdateSolPP2 by
/// Haopeng Wang.
///
/// @author Rayan Dhib

class PhysicalityCFLUpdateSolCorona : public PhysicalityCFLUpdateSol {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit PhysicalityCFLUpdateSolCorona(const std::string& name);

  /// Destructor.
  virtual ~PhysicalityCFLUpdateSolCorona();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

protected:

  /// Largest field strength over the states before they are updated, which
  /// sets the Alfven speed cap. Also runs the one time layout check.
  virtual void beforeUpdate();

  /// Clip the updated states.
  virtual void afterUpdate();

private:

  /// Check that the state ordering this command assumes is the one the
  /// update variable set actually provides.
  void checkStateLayout();

private:

  /// reference density of the corona inner boundary, sets the Alfven speed cap
  CFreal m_rhoBC;

  /// multiplier of the reference Alfven speed in the cap
  CFreal m_alfvenFactor;

  /// temperature floor [K]
  CFreal m_minT;

  /// radius below which the velocity damping applies
  CFreal m_dampingRadius;

  /// mean molecular weight used in the temperature
  CFreal m_mu;

  /// corona normalization, [Bref, rhoref, pref, vref] in SI
  std::vector<CFreal> m_refValues;

  /// largest field strength in the domain before the current update
  CFreal m_bMagMax;

  /// the state layout has been checked against the update variables
  bool m_layoutChecked;

}; // class PhysicalityCFLUpdateSolCorona

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_PhysicalityCFLUpdateSolCorona_hh
