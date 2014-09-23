#ifndef COOLFluiD_Numerics_FluctSplit_NullJacobianFix_hh
#define COOLFluiD_Numerics_FluctSplit_NullJacobianFix_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeJacobianFix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {
      class FluctuationSplitData;

//////////////////////////////////////////////////////////////////////////////

/// This class computes a null jacobian fix to cure entropy violation or carbuncle
/// for compressible flows
/// @author Andrea Lani
class FluctSplit_API NullJacobianFix : public ComputeJacobianFix {
public:

  /// Default constructor without arguments
  NullJacobianFix(const std::string& name);

  /// Default destructor
  virtual ~NullJacobianFix();

  /// Set up
  virtual void setup();

  /// Compute the jacobian fix
  virtual void computeFix(const InwardNormalsData& normalsData,
                          RealVector& delta);

}; // end of class NullJacobianFix

//////////////////////////////////////////////////////////////////////////////

    }  // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NullJacobianFix_hh
