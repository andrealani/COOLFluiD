#ifndef COOLFluiD_Numerics_FluctSplit_PicardJacobStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_PicardJacobStrategy_hh

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
class FluctSplit_API PicardJacobStrategy : public ComputeJacobStrategy {
public:

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer VectorTransformer;

  /// Constructor.
  PicardJacobStrategy(const std::string& name);

  /// Destructor.
  ~PicardJacobStrategy();

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
			   const std::vector<CFuint>&  equationIDs);
  
protected:

  /// Set up private data and data
  void setup();

private:

  /// the single splitter
  Common::SafePtr<Splitter> _splitter;

  /// the system splitter
  Common::SafePtr<Splitter> _sysSplitter;

  /// the scalar splitter
  Common::SafePtr<Splitter> _scalarSplitter;

  // matrix in which storing the temporary value of the jacobian
  // contribution
  RealMatrix _jacob;

  // matrix in which storing the temporary value of the jacobian
  // contribution
  RealMatrix _temp0;

  // matrix in which storing the temporary value of the jacobian
  // contribution
  RealMatrix _temp1;

  // matrix in which storing the temporary value of the jacobian
  // contribution
  RealMatrix _temp2;

  // vector of matrices (blocks) storing the block contributions
  // to the jacobians for each state of the current cell
  // the storage is ROW ORIENTED
  std::vector<RealMatrix*> _jacobBlocks;

  // vector of matrices (blocks) storing the block contributions
  // to the jacobians of the diffusion flux for each state of the current cell
  // the storage is ROW ORIENTED
  std::vector<RealMatrix*> _diffjacob;

}; // class PicardJacobStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PicardJacobStrategy_hh
