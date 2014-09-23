#ifndef COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxCompactApproach_hh
#define COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxCompactApproach_hh

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
 * based on a compact approach
 *
 * @author Kris Van den Abeele
 */
class FaceDiffusiveFluxCompactApproach : public FaceDiffusiveFlux {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,FaceDiffusiveFluxCompactApproach > PROVIDER;

public:  // methods

  /// Constructor
  FaceDiffusiveFluxCompactApproach(const std::string& name) : FaceDiffusiveFlux(name) {}

  /// Destructor
  ~FaceDiffusiveFluxCompactApproach() {}

  /// Compute averaged physical variable in a series of flux points, from left and right states and face normal
  virtual std::vector< RealVector >& computeAvgGradVar(const CFuint iVar,
                                                       std::vector< RealVector* >& lStates,
                                                       std::vector< RealVector* >& rStates,
                                                       const CFuint nbrFlxPnts) = 0;

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FaceDiffusiveFluxCompactApproach";
  }
}; // class FaceDiffusiveFluxCompactApproach

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_FaceDiffusiveFluxCompactApproach_hh

