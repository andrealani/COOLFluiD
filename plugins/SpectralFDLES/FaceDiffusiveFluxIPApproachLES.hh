#ifndef COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxIPApproachLES_hh
#define COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxIPApproachLES_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NSFaceDiffusiveFluxIPApproach.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the averaged diffusive flux on a face,
 * based on the `Interior Penalty' approach
 *
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class FaceDiffusiveFluxIPApproachLES : public NSFaceDiffusiveFluxIPApproach {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,FaceDiffusiveFluxIPApproachLES > PROVIDER;

public:  // methods

  /// Constructor
  FaceDiffusiveFluxIPApproachLES(const std::string& name);

  /// Destructor
  ~FaceDiffusiveFluxIPApproachLES();

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
    return "FaceDiffusiveFluxIPApproachLES";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Set filter width volumes
  void setFilterWidthVolumes(Common::SafePtr< RealVector > filterWidthVolumes)
  {
    m_filterWidthVolumes = filterWidthVolumes;
  }

protected: // data

  /// LES variable set
  Common::SafePtr< LES::LESVarSet > m_lesVarSet;

  /// variable for multiple flux point filter widths on the face
  Common::SafePtr< RealVector > m_filterWidthVolumes;

}; // class FaceDiffusiveFluxIPApproachLES

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxIPApproachLES_hh

