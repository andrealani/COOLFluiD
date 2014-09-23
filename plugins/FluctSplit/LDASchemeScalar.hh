#ifndef COOLFluiD_Numerics_FluctSplit_LDASchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_LDASchemeScalar_hh

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
class FluctSplitScalar_API LDASchemeScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  LDASchemeScalar(const std::string& name);

  /// Default destructor
  ~LDASchemeScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  void distribute(std::vector<RealVector>& residual);

  /// Distribute part of the residual
  void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealVector _phi;

  RealVector _sumKplus;

  RealVector _uTemp;

  RealVector _temp;

}; // end of class LDASchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeScalar_hh
