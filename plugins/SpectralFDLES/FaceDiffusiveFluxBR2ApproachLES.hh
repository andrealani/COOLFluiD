#ifndef COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxBR2ApproachLES_hh
#define COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxBR2ApproachLES_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NSFaceDiffusiveFluxBR2Approach.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the averaged diffusive flux on a face,
 * based on the Bassi-Rebay II approach
 *
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class FaceDiffusiveFluxBR2ApproachLES : public NSFaceDiffusiveFluxBR2Approach {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,FaceDiffusiveFluxBR2ApproachLES > PROVIDER;

public:  // methods

  /// Constructor
  FaceDiffusiveFluxBR2ApproachLES(const std::string& name);

  /// Destructor
  ~FaceDiffusiveFluxBR2ApproachLES();

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
    return "FaceDiffusiveFluxBR2ApproachLES";
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

}; // class FaceDiffusiveFluxBR2ApproachLES

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxBR2ApproachLES_hh

