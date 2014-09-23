#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerNumJacobStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerNumJacobStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeJacobStrategy.hh"
#include "Framework/BlockAccumulator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a fluctuation splitting strategy
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerNumJacobStrategy : public ComputeJacobStrategy {
public:

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Constructor.
   */
  TwoLayerNumJacobStrategy(const std::string& name);

  /**
   * Destructor.
   */
  ~TwoLayerNumJacobStrategy();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeJacobStrategy::configure(args);
  }

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeJacobStrategy::needsSockets();

    result.push_back(&socket_interStates);

    return result;
  }

  /**
   * Add compute the term to add in the jacobian
   */
  void computeJacobianTerm (Framework::GeometricEntity *const cell,
			    const std::vector<RealVector>& residual,
			    Framework::BlockAccumulator *const acc,
			    const std::vector<CFuint>&  equationIDs);
  
protected:

  /**
   * Set up private data and data
   */
  void setup();

  /**
   * Set up private data and data
   */
  void cleanOtherResidual()
  {
    const CFuint n = _otherResidual.size();
    for (CFuint i = 0; i < n; i++) {
      _otherResidual[i] = 0.;
    }
  }

 private:

  // handle to states
  Framework::DataSocketSink<Framework::State*> socket_interStates;

  /// storage for the temporary node residuals
  RealVector _tempRes;

  /// storage for the temporary perturbed states
  std::vector<RealVector> _otherResidual;

  /// storage for the temporary perturbed states
  std::vector<RealVector> _tResidual;

}; // class TwoLayerNumJacobStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerNumJacobStrategy_hh
