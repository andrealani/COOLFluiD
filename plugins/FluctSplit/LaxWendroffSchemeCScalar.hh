#ifndef COOLFluiD_Numerics_FluctSplit_LaxWendroffSchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_LaxWendroffSchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the LaxWendroff scheme for RDS space discretization
/// based on the CRD approach.
/// The implementation was taken as-is from THOR code.
/// @author Tiago Quintino
class FluctSplitScalar_API LaxWendroffSchemeCScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  LaxWendroffSchemeCScalar(const std::string& name);

  /// Default destructor
  ~LaxWendroffSchemeCScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Distribute the residual
  virtual void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private: // data

  RealVector m_sumKabs;

}; // end of class LaxWendroffSchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LaxWendroffSchemeCScalar_hh
