#ifndef COOLFluiD_FluxReconstructionMethod_PositivityIDPMHD_hh
#define COOLFluiD_FluxReconstructionMethod_PositivityIDPMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/BasePositivityIDP.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
      class MHD2DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Invariant-domain-preserving limiter for ideal MHD with GLM divergence
 * cleaning.
 *
 * Conservative state: (rho, rhoU, rhoV, rhoW, Bx, By, Bz, rhoE, phi)
 * Update state (Prim): (rho, u, v, w, Bx, By, Bz, p, phi)
 *
 * The layout is dimension independent: the 2D projection model keeps the z
 * components of velocity and B even though they are typically zero, so both
 * 2D and 3D carry 9 equations. Indices must therefore never be derived from
 * the dimension.
 *
 * The admissible set is {rho > 0, p > 0}. Pressure is concave in the
 * conservative state on {rho > 0}, since it is a linear term minus the
 * quadratic-over-linear perspective function |m|^2/(2 rho) minus the convex
 * |B|^2/2. So the set is convex and the Zhang-Shu scaling applies, with the
 * concavity chord giving a closed-form sufficient theta.
 *
 * @author Rayan Dhib
 */

class PositivityIDPMHD : public BasePositivityIDP {
public:

  explicit PositivityIDPMHD(const std::string& name);

  virtual ~PositivityIDPMHD();

  static void defineConfigOptions(Config::OptionList& options);

protected: // physics hooks

  virtual void constraintsAtPoint(const RealVector& cons,
                                  CFreal& rho, CFreal& p, CFreal& B2) const;

  virtual void consToUpdate(const RealVector& cons, RealVector& update) const;

  virtual const std::vector< CFuint >& scaledIndices(ScaleMode mode) const;

  virtual void setupPhysics();

protected: // data

  /// gamma - 1
  CFreal m_gammaMinusOne;

  /// conservative indices scaled in Hydro mode: rho, momentum, energy
  std::vector< CFuint > m_hydroIndices;

  /// conservative indices scaled in Full mode: everything
  std::vector< CFuint > m_fullIndices;

  /// index of the energy component in the conservative state
  CFuint m_iE;

  /// index of the first magnetic field component
  CFuint m_iB;

}; // class PositivityIDPMHD

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PositivityIDPMHD_hh
