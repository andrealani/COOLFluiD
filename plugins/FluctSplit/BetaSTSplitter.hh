#ifndef COOLFluiD_Numerics_FluctSplit_BetaSTSplitter_hh
#define COOLFluiD_Numerics_FluctSplit_BetaSTSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "SourceTermSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a splitter for source term based on lumped Galerkin
/// @author Andrea Lani
class FluctSplit_API BetaSTSplitter : public SourceTermSplitter {
public:

  /// Default constructor without arguments
  BetaSTSplitter(const std::string& name);

  /// Default destructor
  virtual ~BetaSTSplitter();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Compute the source term
  virtual void computeSourceTerm(const InwardNormalsData& normalsData);

private:

  /// Distribute the residual
  virtual void distributeExceptBStates(std::vector<RealVector>& residual);

private:

  /// temporary source term
  RealVector _tmpSource;

  /// flag telling if to exclude the distribution to boundary states
  bool _excludeBStates;

}; // end of class BetaSTSplitter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BetaSTSplitter_hh
