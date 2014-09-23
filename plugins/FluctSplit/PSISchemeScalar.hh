#ifndef COOLFluiD_Numerics_FluctSplit_PSISchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_PSISchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents the N scheme for RDS space discretization
  /// @author Andrea Lani
  /// @author Tiago Quintino
class FluctSplitScalar_API PSISchemeScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  PSISchemeScalar(const std::string& name);

  /// Default destructor
  ~PSISchemeScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  void distribute(std::vector<RealVector>& residual);

  /// Distribute part of the residual
  void distributePart(std::vector<RealVector>& residual);

private:

  RealVector _sumKminU;

  RealVector _sumKmin;

  RealVector _uTemp;

  RealVector _uDiff;

  RealVector _phy;

  RealVector _sumBeta;

  RealVector _invCoeff;

  RealVector _temp;

}; // end of class PSISchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PSISchemeScalar_hh
