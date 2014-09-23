#ifndef COOLFluiD_Numerics_FluctSplit_NullComputeDiffusiveTerm_hh
#define COOLFluiD_Numerics_FluctSplit_NullComputeDiffusiveTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class computes the diffusive flux corresponding to the Navier Stokes
/// physical model
/// @author Andrea Lani
class FluctSplit_API NullComputeDiffusiveTerm : public ComputeDiffusiveTerm {
public:

  /// Constructor
  NullComputeDiffusiveTerm(const std::string& name);

  /// Default destructor
  ~NullComputeDiffusiveTerm();

  /// Set up private data to prepare the simulation
  void setup();

  /// Set the update variable set
  void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar);
  /// Set the update variable set
  void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet>
  	       updateVar);

  /// Compute the diffusive term flux in the current cell
  void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
      std::vector<RealVector>& result,
      const bool updateCoeff);

  /// Compute the diffusive term flux in the current cell
  bool isNull() const
  {
    return true;
  }

  void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob){
    throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalarNullDifterm::computePicardJacob()");
  }
}; // end of class NullComputeDiffusiveTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NullComputeDiffusiveTerm_hh
