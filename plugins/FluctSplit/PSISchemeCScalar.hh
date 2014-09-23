#ifndef COOLFluiD_Numerics_FluctSplit_PSISchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_PSISchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents the PSI scheme for CRD space discretization
  /// @author Andrea Lani
class FluctSplitScalar_API PSISchemeCScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  PSISchemeCScalar(const std::string& name);

  /// Default destructor
  ~PSISchemeCScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Distribute part of the residual
  virtual void distributePart(std::vector<RealVector>& residual);

private:

  RealVector _sumKplus;
  
  RealVector _sumKplusU;

  RealVector _uInflow;

  RealVector _uDiff;

  RealVector _sumBeta;

  RealVector _invCoeff;

  RealVector _temp;

}; // end of class PSISchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PSISchemeCScalar_hh
