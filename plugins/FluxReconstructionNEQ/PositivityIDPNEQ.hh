// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_PositivityIDPNEQ_hh
#define COOLFluiD_FluxReconstructionMethod_PositivityIDPNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/BasePositivityIDP.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * IDP (Zhang-Shu) positivity limiter for NEQ, keeping every partial density
 * positive.
 *
 * Scope: species densities only. The thermal constraint is deliberately left
 * out, see below. Negative partial densities are the failure mode that breaks
 * the chemistry library, so this removes the crash without weakening
 * conservation.
 *
 * Why species only, and why that costs nothing:
 * the first nbSpecies components of the update set (RhoivtTv) and of the
 * conservative set are the same partial densities. Scaling only those means
 * the limited components need no variable transformation at all, so
 * writeBackStates() is overridden to write them straight back and the
 * cons -> update mapping is never needed. Species mass stays exactly conserved,
 * since the scaling is affine about the metric-weighted mean of a conservative
 * quantity. Momentum and energy are not exactly conserved: u, v, T and Tv are
 * held fixed while rho changes, so rho*u and rho*E move slightly, and only in
 * the cells where the limiter fires.
 *
 * Why no thermal constraint: the quantity that is positive iff T > 0 is the
 * SENSIBLE internal energy, not the total. For a N/N2 mixture the formation
 * energies are large and positive, so total internal energy stays positive at
 * T = 0 and a constraint built on it would be silently inert. Sensible energy
 * needs species formation energies, which would come from
 * PhysicalChemicalLibrary::ExtraData, and PlatoLibrary fills no ExtraData at
 * all. The pressure constraint is therefore forced off in setupPhysics().
 *
 * Adding it later means writing the cons -> RhoivtTv inversion (a Newton solve
 * on T and Tv against the library) and constraining sensible energy.
 *
 * @author Rayan Dhib
 */
class PositivityIDPNEQ : public BasePositivityIDP {
public:

  explicit PositivityIDPNEQ(const std::string& name);

  virtual ~PositivityIDPNEQ();

  static void defineConfigOptions(Config::OptionList& options);

protected: // physics hooks

  /// rho is the smallest partial density; p is not constrained
  virtual void constraintsAtPoint(const RealVector& cons,
                                  CFreal& rho, CFreal& p, CFreal& B2) const;

  /// Not supported, and not reachable: writeBackStates() is overridden so the
  /// general cons -> update mapping is never needed. Throws if ever called.
  virtual void consToUpdate(const RealVector& cons, RealVector& update) const;

  /// the species density components, in both scaling modes
  virtual const std::vector< CFuint >& scaledIndices(ScaleMode mode) const;

  virtual void setupPhysics();

protected: // machinery

  /// Write back only the species densities, which are the same entries in the
  /// conservative and in the update set
  virtual void writeBackStates();

protected: // data

  /// number of species
  CFuint m_nbSpecies;

  /// indices of the partial densities, the only scaled components
  std::vector< CFuint > m_speciesIndices;

}; // class PositivityIDPNEQ

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PositivityIDPNEQ_hh
