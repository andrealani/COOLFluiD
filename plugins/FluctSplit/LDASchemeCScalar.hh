#ifndef COOLFluiD_Numerics_FluctSplit_LDASchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_LDASchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the LDA scheme for RDS space discretization
/// based on the CRD approach
/// @author Tiago Quintino
class FluctSplitScalar_API LDASchemeCScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  LDASchemeCScalar(const std::string& name);

  /// Default destructor
  ~LDASchemeCScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Distribute the residual
  virtual void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealVector m_sumKplus;

  RealVector m_uTemp;

}; // end of class LDASchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeCScalar_hh
