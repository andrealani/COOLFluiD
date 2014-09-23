#ifndef COOLFluiD_Numerics_FluctSplit_NSchemeScalarT_hh
#define COOLFluiD_Numerics_FluctSplit_NSchemeScalarT_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MacroLET.hh"
#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the N scheme for RDS space discretization
/// @author Andrea Lani
/// @author Tiago Quintino
template <CFuint N>
class FluctSplitScalar_API NSchemeScalarT : public RDS_SplitterScalar {
public:

  typedef NSchemeScalarT<N> SELF;

  /// Default constructor.
  NSchemeScalarT(const std::string& name);

  /// Default destructor
  ~NSchemeScalarT();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  void distribute(std::vector<RealVector>& residual);

private:

  HPVECTOR(N) _sumKminU;

  HPVECTOR(N) _sumKmin;

  HPVECTOR(N) _uTemp;

  HPVECTOR(N) _uMin;

  HPVECTOR(N) _temp;

  HPVECTORDYN _tStateTmp;

  HPVECTORDYN _kTmp;

}; // end of class NSchemeScalarT

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NSchemeScalarT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeScalarT_hh
