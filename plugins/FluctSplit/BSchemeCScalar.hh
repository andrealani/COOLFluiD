#ifndef COOLFluiD_Numerics_FluctSplit_BSchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_BSchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/BSchemeBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the Blended scheme for RDS space discretization
/// based on the CRD approach
/// @author Tiago Quintino
class FluctSplitScalar_API BSchemeCScalar : public BSchemeBase<RDS_SplitterScalar> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit BSchemeCScalar(const std::string& name);

  /// Virtual destructor
  virtual ~BSchemeCScalar();

  /// Configure this object with user defined parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup this object with data depending on the mesh
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  ///  Distribute part of the residual
  virtual void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

protected: // functions

  /// Compute the blending coefficients
  virtual void computeBlendingCoeff();

  /**
   * Add an isotropic dissipative term.
   */
  virtual void addExtraDissipation();

private:

  RealVector m_sumKplusU;

  RealVector m_sumKplus;

  RealVector m_phiLDA;

  RealVector m_uInFlow;

  RealVector m_temp;

  RealVector m_uMin;

}; // end of class BSchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BSchemeCScalar_hh
