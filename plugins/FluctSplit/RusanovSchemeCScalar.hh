#ifndef COOLFluiD_Numerics_FluctSplit_RusanovSchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_RusanovSchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the Rusanov scheme for RDS space discretization
/// based on the CRD approach
/// @author Tiago Quintino
class FluctSplitScalar_API RusanovSchemeCScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  RusanovSchemeCScalar(const std::string& name);

  /// Default destructor
  ~RusanovSchemeCScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Distribute the residual
  virtual void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private: // data

  // diffusion coefficient
  RealVector m_alpha;

  RealVector m_sumUmin;

}; // end of class RusanovSchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RusanovSchemeCScalar_hh
