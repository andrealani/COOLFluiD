#ifndef COOLFluiD_FluxReconstructionMethod_HLLFlux_hh
#define COOLFluiD_FluxReconstructionMethod_HLLFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"
#include <stdio.h>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent an HLL interface flux
 *
 * @author Ray Vandenhoeck
 */
class HLLFlux : public RiemannFlux {

public:  // methods

  /// Constructor
  HLLFlux(const std::string& name);

  /// Destructor
  ~HLLFlux();
 
  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState,
				  Framework::State& rState,
				  const RealVector& normal);

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);


  /// Gets the Class name
  static std::string getClassName()
  {
    return "HLLFlux";
  }

  /// Set up private data and data
  virtual void setup();
  
  /// Unset up private data and data
  virtual void unsetup();

private: // data

  /// array storing the sum of the right flux
  RealVector   m_rightFlux;

 /// array storing the sum of the left flux
  RealVector   m_leftFlux;

  /// array storing the temporary right eigenvalues
  RealVector    m_rightEv;

  /// array storing the temporary left eigenvalues
  RealVector    m_leftEv;

  /// temporary unit normal
  RealVector    m_tempUnitNormal;

  /// pre-allocated left solution state (avoids heap allocation per flux call)
  RealVector    m_lSolState;

  /// pre-allocated right solution state (avoids heap allocation per flux call)
  RealVector    m_rSolState;


}; // class HLLFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_HLLFlux_hh
