#ifndef COOLFluiD_Numerics_FluctSplit_BSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_BSchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/BSchemeBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents the B scheme for RDS space discretization
  /// @author Andrea Lani
  /// @author Tiago Quintino
class FluctSplitScalar_API BSchemeScalar : public BSchemeBase<RDS_SplitterScalar> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit BSchemeScalar(const std::string& name);

  /// Destructor
  virtual ~BSchemeScalar();

  /// Configure this object with user defined parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup this object with data depending on the mesh
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  ///  Distribute part of the residual
  virtual void distributePart(std::vector<RealVector>& residual);

protected: // functions

  /// Compute the blending coefficients
  virtual void computeBlendingCoeff();

  /**
   * Add an isotropic dissipative term.
   */
  virtual void addExtraDissipation();

private:

  RealVector _sumKminU;

  RealVector _sumKmin;

  RealVector _sumKplus;

  RealVector _uInflow;

  RealVector _uDiff;

  RealVector _phi;

  RealVector _phiLDA;

  RealVector _temp;

}; // end of class BSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BSchemeScalar_hh
