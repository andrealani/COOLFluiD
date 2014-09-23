#ifndef COOLFluiD_Numerics_FluctSplit_LumpedGalerkinSTSplitter_hh
#define COOLFluiD_Numerics_FluctSplit_LumpedGalerkinSTSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "SourceTermSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a splitter for source term based on lumped Galerkin
/// @author Andrea Lani
class FluctSplit_API LumpedGalerkinSTSplitter : public SourceTermSplitter {
public:

  /// Default constructor without arguments
  LumpedGalerkinSTSplitter(const std::string& name);

  /// Default destructor
  virtual ~LumpedGalerkinSTSplitter();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Compute the source term
  virtual void computeSourceTerm(const InwardNormalsData& normalsData);

}; // end of class LumpedGalerkinSTSplitter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LumpedGalerkinSTSplitter_hh
