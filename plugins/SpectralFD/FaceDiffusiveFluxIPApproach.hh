#ifndef COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxIPApproach_hh
#define COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxIPApproach_hh

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
 * based on the `Interior Penalty' approach
 *
 * @author Kris Van den Abeele
 */
class FaceDiffusiveFluxIPApproach : public FaceDiffusiveFluxCompactApproach {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,FaceDiffusiveFluxIPApproach > PROVIDER;

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  FaceDiffusiveFluxIPApproach(const std::string& name);

  /// Destructor
  ~FaceDiffusiveFluxIPApproach();

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
    return "FaceDiffusiveFluxIPApproach";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

protected: // data

  /// parameter of the IP approach determining the amount of damping
  CFreal m_alpha;

  /// solution polynomial order
  CFreal m_solPolyOrderFac;

  /// averaged solution
  RealVector m_avgSol;

  /// averaged gradients
  std::vector< RealVector* > m_avgGrads;

  /// left damping terms
//   std::vector< RealVector > m_dampTermL;

  /// right damping terms
//   std::vector< RealVector > m_dampTermR;

}; // class FaceDiffusiveFluxIPApproach

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxIPApproach_hh

