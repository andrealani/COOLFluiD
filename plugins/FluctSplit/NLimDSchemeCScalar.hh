#ifndef COOLFluiD_Numerics_FluctSplit_NLimDSchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_NLimDSchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the N Limited scheme for RDS space discretization
/// based on the CRD approach but with the addition of
/// an artificial dissipation based on Abgrall's paper on JCP 2005
/// "Essential non-oscillatory Residual Distribution Schemes for hyperbolic problems."
/// @author Tiago Quintino
class FluctSplitScalar_API NLimDSchemeCScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  NLimDSchemeCScalar(const std::string& name);

  /// Default destructor
  ~NLimDSchemeCScalar();

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

  /// the residual of the standard N scheme before limiting is applied
  std::vector<RealVector> m_phiN;

  /// the beta coefficient relative to the unlimited N scheme
  std::vector<RealVector> m_betaN;

  /// the positive values of the beta coefficients if the unlimited N scheme
  std::vector<RealVector> m_betaNPlus;

  /// the sum of the m_betaPlus over all the states
  RealVector m_sumBetaNPlus;

  /// the final distribution coefficient for the limited N scheme
  std::vector<RealVector> m_beta;

  /// the theta function for detecting the presence of discontinuities
  RealVector m_theta;

  /// the average state in the cell
  RealVector m_ubar;

}; // end of class NLimDSchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NLimDSchemeCScalar_hh
