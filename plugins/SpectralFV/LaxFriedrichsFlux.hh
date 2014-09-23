#ifndef COOLFluiD_Numerics_SpectralFV_LaxFriedrichsFlux_hh
#define COOLFluiD_Numerics_SpectralFV_LaxFriedrichsFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFV/SpectralFVMethodData.hh"
#include "SpectralFV/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a Lax-Friedrichs/Rusanov flux
 *
 * @author Kris Van den Abeele
 */
class LaxFriedrichsFlux : public RiemannFlux {

public:  // methods

  /// Constructor
  LaxFriedrichsFlux(const std::string& name);

  /// Destructor
  ~LaxFriedrichsFlux();

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);

  /// Compute the Riemann flux, from a vector of left and right states and normal vectors
  /// States data must be computed before calling this function!!!
  virtual std::vector< RealVector >& computeFlux(std::vector< Framework::State* >& lState,
                                                 std::vector< Framework::State* >& rState,
                                                 const std::vector< RealVector >& normal,
                                                 const CFuint nbrFlxPnts);

  /// Compute the numerical damping term in the Riemann flux,
  /// from left and right states and a normal vector
  virtual RealVector& computeNumDamping(Framework::State& lState, RealVector& lExtraVars,
                                        Framework::State& rState, RealVector& rExtraVars,
                                        const RealVector& normal);

  /// Compute the numerical damping term in the Riemann flux,
  /// from a vector of left and right states and normal vectors
  /// States data must be computed before calling this function!!!
  virtual std::vector< RealVector >& computeNumDamping(std::vector< Framework::State* >& lState,
                                                       std::vector< Framework::State* >& rState,
                                                       const std::vector< RealVector >& normal,
                                                       const CFuint nbrFlxPnts);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LaxFriedrichsFlux";
  }

  /// Set up private data and data
  virtual void setup();

private: // data

  /// variable set transformer from update to solution variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_updateToSolutionVarTrans;

  /// array storing the sum of the right and left flux
  RealVector  m_sumFlux;

  /// vector storing the left and right update states
  std::vector<Framework::State*> m_updateStates;

  /// vector storing the left and right update states in one flux point
  std::vector<Framework::State*> m_updateStatesInFlxPnt;

  /// vector storing left and right extra variables in one flux point
  std::vector< RealVector* > m_extraVarsInFlxPnt;

  /// vector storing the left and right solution states
  std::vector<Framework::State*>* m_solStates;

  /// vector storing the left and right solution states in one flux point
  std::vector<Framework::State*> m_solStatesInFlxPnt;

}; // class LaxFriedrichsFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_LaxFriedrichsFlux_hh
