#ifndef COOLFluiD_Numerics_SpectralFV_RiemannFlux_hh
#define COOLFluiD_Numerics_SpectralFV_RiemannFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Riemann solver
 *
 * @author Kris Van den Abeele
 * @author Tiago Quintino
 */
class RiemannFlux : public SpectralFVMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFVMethodData,RiemannFlux > PROVIDER;

public:  // methods

  /// Constructor
  RiemannFlux(const std::string& name);

  /// Destructor
  ~RiemannFlux();

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal) = 0;

  /// Compute the Riemann flux, from a vector of left and right states and normal vectors
  /// States data must be computed before calling this function!!!
  virtual std::vector< RealVector >& computeFlux(std::vector< Framework::State* >& lState,
                                                 std::vector< Framework::State* >& rState,
                                                 const std::vector< RealVector >& normal,
                                                 const CFuint nbrFlxPnts) = 0;

  /// Compute the numerical damping term in the Riemann flux,
  /// from left and right states and a normal vector
  virtual RealVector& computeNumDamping(Framework::State& lState, RealVector& lExtraVars,
                                        Framework::State& rState, RealVector& rExtraVars,
                                        const RealVector& normal) = 0;

  /// Compute the numerical damping term in the Riemann flux,
  /// from a vector of left and right states and normal vectors
  /// States data must be computed before calling this function!!!
  virtual std::vector< RealVector >& computeNumDamping(std::vector< Framework::State* >& lState,
                                                       std::vector< Framework::State* >& rState,
                                                       const std::vector< RealVector >& normal,
                                                       const CFuint nbrFlxPnts) = 0;

  /// Gets the Class name
  static std::string getClassName()
  {
    return "RiemannFlux";
  }
 
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
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

  /// variable for numerical damping
  RealVector m_numDamping;

  /// variable for multiple Riemann fluxes
  std::vector< RealVector > m_multiRFlux;

  /// variable for multiple numerical dampings
  std::vector< RealVector > m_multiNumDamping;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class RiemannFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_RiemannFlux_hh

