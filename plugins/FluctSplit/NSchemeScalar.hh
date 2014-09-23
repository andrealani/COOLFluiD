#ifndef COOLFluiD_Numerics_FluctSplit_NSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_NSchemeScalar_hh

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
class FluctSplitScalar_API NSchemeScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  NSchemeScalar(const std::string& name);

  /// Default destructor
  ~NSchemeScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  void distribute(std::vector<RealVector>& residual);

  /// Distribute part of the residual
  void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealVector _sumKminU;

  RealVector _sumKmin;

  RealVector _uTemp;

  RealVector _uMin;

  RealVector _temp;

}; // end of class NSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeScalar_hh
