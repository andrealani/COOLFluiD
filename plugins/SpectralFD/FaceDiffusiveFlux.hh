#ifndef COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFlux_hh
#define COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/DiffusiveVarSet.hh"

#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the averaged diffusive flux on a face
 *
 * @author Kris Van den Abeele
 */
class FaceDiffusiveFlux : public SpectralFDMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,FaceDiffusiveFlux > PROVIDER;

public:  // methods

  /// Constructor
  FaceDiffusiveFlux(const std::string& name);

  /// Destructor
  ~FaceDiffusiveFlux();

  /// Compute averaged gradient variables in a series of flux points, from left and right states and face normal
  virtual std::vector< RealVector >& computeAvgGradVars(std::vector< RealVector* >& lStates,
                                                        std::vector< RealVector* >& rStates,
                                                        const CFuint nbrFlxPnts) = 0;

  /// Compute the averaged diffusive flux in a series of flux points,
  /// from left and right gradients and states and a normal vector
  virtual std::vector< RealVector >& computeDiffFlux(std::vector< std::vector< RealVector* >* >& lGrads,
                                                     std::vector< std::vector< RealVector* >* >& rGrads,
                                                     std::vector< RealVector* >& lStates,
                                                     std::vector< RealVector* >& rStates,
                                                     const std::vector< CFreal >& faceInvCharLength,
                                                     const std::vector< RealVector >& normal,
                                                     const CFuint nbrFlxPnts) = 0;

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FaceDiffusiveFlux";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Set the maximum number of points the Riemann flux will be evaluated in simultaneously
  void setMaxNbrFlxPnts(const CFuint maxNbrFlxPnts)
  {
    m_maxNbrFlxPnts = maxNbrFlxPnts;
  }

protected: // data

  /// diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;

  /// maximum number of flux points the diffusive face flux has to be evaluated for simultaneously
  CFuint m_maxNbrFlxPnts;

  /// variable for multiple diffusive face fluxes
  std::vector< RealVector > m_multiDiffFlux;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class FaceDiffusiveFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFlux_hh

