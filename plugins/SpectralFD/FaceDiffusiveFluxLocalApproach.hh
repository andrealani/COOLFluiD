#ifndef COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxLocalApproach_hh
#define COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxLocalApproach_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/FaceDiffusiveFlux.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the averaged diffusive flux on a face,
 * based on the `local DG' approach
 *
 * @author Kris Van den Abeele
 */
class FaceDiffusiveFluxLocalApproach : public FaceDiffusiveFlux {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,FaceDiffusiveFluxLocalApproach > PROVIDER;

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  FaceDiffusiveFluxLocalApproach(const std::string& name);

  /// Destructor
  ~FaceDiffusiveFluxLocalApproach();

  /// Compute averaged gradient variables in a series of flux points, from left and right states and face normal
  virtual std::vector< RealVector >& computeAvgGradVars(std::vector< RealVector* >& lStates,
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
    return "FaceDiffusiveFluxLocalApproach";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

protected: // data

  /// parameter of the local approach determining the amount of damping
  CFreal m_alpha;

  /// parameter of the local approach determining the wheight of left and right states
  CFreal m_beta;

  /// parameter of the local approach determining the wheight of left and right states
  CFreal m_oEminusBeta;

  /// averaged solution
  RealVector m_avgSol;

  /// averaged gradients
  std::vector< RealVector* > m_avgGrads;

}; // class FaceDiffusiveFluxLocalApproach

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxLocalApproach_hh

