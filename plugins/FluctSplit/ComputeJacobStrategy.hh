#ifndef COOLFluiD_Numerics_FluctSplit_ComputeJacobStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeJacobStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Framework/DataSocketSink.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a strategy to compute the global jacobian matrix
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API ComputeJacobStrategy : public Framework::MethodStrategy<FluctuationSplitData> {

public: // typedefs

  typedef Framework::BaseMethodStrategyProvider<FluctuationSplitData,ComputeJacobStrategy> PROVIDER;
  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer       VectorTransformer;

  /// Constructor.
  ComputeJacobStrategy(const std::string& name);

  /// Destructor.
  virtual ~ComputeJacobStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Framework::MethodStrategy<FluctuationSplitData>::configure(args);
  }

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;

    return result;
  }

  /// Add compute the term to add in the jacobian
  virtual void computeJacobianTerm (Framework::GeometricEntity *const cell,
				    const std::vector<RealVector>& residual,
				    Framework::BlockAccumulator *const acc,
				    const std::vector<CFuint>& equationIDs) = 0;
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeJacobStrategy";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /// Set up private data and data
  virtual void setup();

protected: // data

  /// acquaintance of the FluctuationSplitStrategy
  Common::SafePtr<FluctuationSplitStrategy> _fsStrategy;

/// acquaintance of the artificial diffusion strategy Strategy
  Common::SafePtr<ArtificialDiffusionStrategy> _adStrategy;

  /// Transformer from Solution to Distribution Variables
  Common::SafePtr<MatrixTransformer> _solutionToDistMatTrans;

  /// Transformer from Distribution to Solution Variables
  Common::SafePtr<MatrixTransformer> _distToSolutionMatTrans;

  /// Transformer from Linear to Distribution Variables
  Common::SafePtr<MatrixTransformer> _linearToDistMatTrans;

  /// Transformer from Linear to Solution in Solution Variables
  Common::SafePtr<MatrixTransformer> _solutionToLinearInUpdateMatTrans;

  /// Matrix Transformer from Solution to Linear Variables
  Common::SafePtr<MatrixTransformer> _solutionToLinearMatTrans;

  /// Transformer from Update to Linear Variables
  Common::SafePtr<VectorTransformer> _updateToLinearVecTrans;

  /// Matrix Transformer from Update to Solution  in Update Variables
  Common::SafePtr<MatrixTransformer> _updateToSolutionInUpdateMatTrans;

  /// storage for the temporary diffusive residual
  std::vector<RealVector> _diffResidual;

  /// tells if there is a diffusive term to compute
  bool _hasDiffusiveTerm;

}; // class ComputeJacobStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeJacobStrategy_hh
