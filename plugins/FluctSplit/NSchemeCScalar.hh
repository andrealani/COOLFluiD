#ifndef COOLFluiD_Numerics_FluctSplit_NSchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_NSchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the N scheme for RDS space discretization
/// based on the CRD approach
/// @author Andrea Lani
class FluctSplitScalar_API NSchemeCScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  NSchemeCScalar(const std::string& name);

  /// Default destructor
  ~NSchemeCScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Distribute the residual
  virtual void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealVector m_sumKplusU;

  RealVector m_sumKplus;

  RealVector m_uTemp;

  RealVector m_uMin;

  RealVector m_temp;

}; // end of class NSchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeCScalar_hh
