#ifndef COOLFluiD_Numerics_FluctSplit_IncludeInFluxSTSplitter_hh
#define COOLFluiD_Numerics_FluctSplit_IncludeInFluxSTSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "SourceTermSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a splitter for source term based on lumped Galerkin
/// @author Andrea Lani
class FluctSplit_API IncludeInFluxSTSplitter : public SourceTermSplitter {
public:

  /// Default constructor without arguments
  IncludeInFluxSTSplitter(const std::string& name);

  /// Default destructor
  virtual ~IncludeInFluxSTSplitter();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Compute the source term
  virtual void computeSourceTerm(const InwardNormalsData& normalsData);

private:

  /// temporary source term
  RealVector _tmpSource;

}; // end of class IncludeInFluxSTSplitter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_IncludeInFluxSTSplitter_hh
