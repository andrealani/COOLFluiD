#ifndef COOLFluiD_FluxReconstructionMethod_HLLDecEFlux_hh
#define COOLFluiD_FluxReconstructionMethod_HLLDecEFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * HLL Riemann flux with optional decomposed-energy extra diffusion for MHD.
 * FR port of FiniteVolumeMHD/HLLDecEFlux.
 *
 * Adds three switches on top of plain HLL:
 *   UseAlpha          : replace the (amax*amin)/(amax-amin) factor by
 *                       2*alpha*(amax*amin)/(amax-amin),
 *                       alpha = max(|amin|,|amax|)/(amax-amin)
 *   AddLax            : extra plasma-beta-modulated Lax dissipation on the
 *                       energy slot (index 7, MHD3DProjection layout)
 *   2Dornot           : 2D-vs-3D switch. When true (2D case), the plasma-beta
 *                       weighting Q_Betafactor is forced to 1 and the extra
 *                       dissipation is applied unconditionally. When false
 *                       (3D), the dimensional N/m^2 reference values are
 *                       used to form beta and Q_Betafactor = tanh(1/(100*beta)).
 *   Extra_kappaLax_E  : magnitude of the extra dissipation
 *
 * Assumes MHD3DProjection variable ordering: (rho, rhoU, rhoV, rhoW,
 * Bx, By, Bz, rhoE, phi) so that B is at [4,5,6] and rhoE is at [7].
 *
 * @author Rayan Dhib
 * @author Haopeng Wang (FV original)
 */
class HLLDecEFlux : public RiemannFlux {

public:  // methods

  static void defineConfigOptions(Config::OptionList& options);

  HLLDecEFlux(const std::string& name);

  ~HLLDecEFlux();

  virtual RealVector& computeFlux(Framework::State& lState,
                                  Framework::State& rState,
                                  const RealVector& normal);

  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);

  static std::string getClassName()
  {
    return "HLLDecEFlux";
  }

  virtual void setup();

private: // data

  /// right flux
  RealVector m_rightFlux;

  /// left flux
  RealVector m_leftFlux;

  /// right eigenvalues
  RealVector m_rightEv;

  /// left eigenvalues
  RealVector m_leftEv;

  /// left solution state (pre-allocated; transform() returns a shared buffer)
  RealVector m_lSolState;

  /// right solution state
  RealVector m_rSolState;

  /// extra diffusion magnitude on the energy slot
  CFreal m_kappaLax_E;

  /// use modified HLL diffusion coefficient
  bool m_useAlpha;

  /// add extra plasma-beta-modulated Lax dissipation on energy
  bool m_addLax;

  /// 2D-vs-3D case switch; true forces Q_Betafactor = 1
  bool m_2Dornot;

}; // class HLLDecEFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_HLLDecEFlux_hh
