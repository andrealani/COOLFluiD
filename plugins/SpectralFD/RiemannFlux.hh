#ifndef COOLFluiD_Numerics_SpectralFD_RiemannFlux_hh
#define COOLFluiD_Numerics_SpectralFD_RiemannFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Riemann solver
 *
 * @author Kris Van den Abeele
 * @author Tiago Quintino
 */
class RiemannFlux : public SpectralFDMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,RiemannFlux > PROVIDER;

public:  // methods

  /// Constructor
  RiemannFlux(const std::string& name);

  /// Destructor
  ~RiemannFlux();

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal) = 0;
                                  
  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState,
                                  Framework::State& rState,
                                  const RealVector& normal);

  /// Compute the Riemann flux, from a vector of left and right states and normal vectors
  /// States data must be computed before calling this function!!!
  virtual std::vector< RealVector >& computeFlux(std::vector< Framework::State* >& lState,
                                                 std::vector< Framework::State* >& rState,
                                                 const std::vector< RealVector >& normal,
                                                 const CFuint nbrFlxPnts);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "RiemannFlux";
  }

  /// Set up private data and data
  virtual void setup();

  /// Set the maximum number of points the Riemann flux will be evaluated in simultaneously
  void setMaxNbrFlxPnts(const CFuint maxNbrFlxPnts)
  {
    m_maxNbrFlxPnts = maxNbrFlxPnts;
  }
  
protected: // data

  /// maximum number of quadrature (or other) points the Riemann flux has to be evaluated for simultaneously
  CFuint m_maxNbrFlxPnts;

  /// variable for Riemann flux
  RealVector m_rFlux;

  /// variable for multiple Riemann fluxes
  std::vector< RealVector > m_multiRFlux;

  /// number of equations in the physical model
  CFuint m_nbrEqs;
  
  /// temporary vector storing the left and right solution states
  std::vector<Framework::State*> m_solStates;
  
  /// temporary vector storing the left and right physical data
  std::vector<RealVector> m_pData;
  
  /// temporary vector storing the left and right update states
  std::vector<Framework::State*> m_updateStates;

}; // class RiemannFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_RiemannFlux_hh

