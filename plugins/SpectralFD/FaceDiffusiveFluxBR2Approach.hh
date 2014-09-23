#ifndef COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxBR2Approach_hh
#define COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxBR2Approach_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/FaceDiffusiveFluxCompactApproach.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the averaged diffusive flux on a face,
 * based on the Bassi-Rebay II approach
 *
 * @author Kris Van den Abeele
 */
class FaceDiffusiveFluxBR2Approach : public FaceDiffusiveFluxCompactApproach {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,FaceDiffusiveFluxBR2Approach > PROVIDER;

public:  // methods

  /// Constructor
  FaceDiffusiveFluxBR2Approach(const std::string& name);

  /// Destructor
  ~FaceDiffusiveFluxBR2Approach();

  /// Compute averaged states in a series of flux points, from left and right states and face normal
  std::vector< RealVector >& computeAvgGradVars(std::vector< RealVector* >& lStates,
                                                        std::vector< RealVector* >& rStates,
                                                        const CFuint nbrFlxPnts);

  /// Compute averaged physical variable in a series of flux points, from left and right states and face normal
  std::vector< RealVector >& computeAvgGradVar(const CFuint iVar,
                                               std::vector< RealVector* >& lStates,
                                               std::vector< RealVector* >& rStates,
                                               const CFuint nbrFlxPnts);

  /// Compute the averaged diffusive flux in a series of flux points,
  /// from left and right gradients and states and a normal vector
  virtual std::vector< RealVector >& computeDiffFlux(std::vector< std::vector< RealVector* >* >& lGrads,
                                                     std::vector< std::vector< RealVector* >* >& rGrads,
                                                     std::vector< RealVector* >& lStates,
                                                     std::vector< RealVector* >& rStates,
                                                     const std::vector< CFreal >& faceInvCharLength,
                                                     const std::vector< RealVector >& normal,
                                                     const CFuint nbrFlxPnts);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FaceDiffusiveFluxBR2Approach";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

protected: // data

  /// solution polynomial order
  CFreal m_solPolyOrderFac;

  /// averaged solution
  RealVector m_avgSol;

  /// averaged gradients
  std::vector< RealVector* > m_avgGrads;

}; // class FaceDiffusiveFluxBR2Approach

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxBR2Approach_hh

