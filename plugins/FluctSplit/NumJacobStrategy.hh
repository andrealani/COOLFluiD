#ifndef COOLFluiD_Numerics_FluctSplit_NumJacobStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_NumJacobStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeJacobStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Andrea Lani
class FluctSplit_API NumJacobStrategy : public ComputeJacobStrategy {
public:

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer VectorTransformer;

  /// Constructor.
  NumJacobStrategy(const std::string& name);

  /// Destructor.
  ~NumJacobStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeJacobStrategy::configure(args);
  }

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeJacobStrategy::needsSockets();

    return result;
  }

  /// Add compute the term to add in the jacobian
  void computeJacobianTerm(Framework::GeometricEntity *const cell,
			   const std::vector<RealVector>& residual,
			   Framework::BlockAccumulator *const acc,
			   const std::vector<CFuint>& equationIDs);
  
protected:

  /// Set up private data and data
  void setup();

  /// Set up private data and data
  void cleanOtherResidual()
  {
    fill(_otherResidual.begin(),_otherResidual.end(),0.0);
  }

private:

  /// storage for the temporary node residuals
  RealVector _tempRes;

  /// storage for the temporary perturbed states
  std::vector<RealVector> _otherResidual;

  /// storage for the temporary perturbed states
  std::vector<RealVector> _tResidual;

}; // class NumJacobStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NumJacobStrategy_hh
