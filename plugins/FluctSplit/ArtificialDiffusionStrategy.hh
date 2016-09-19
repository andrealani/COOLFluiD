#ifndef COOLFluiD_Numerics_FluctSplit_ArtificialDiffusionStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_ArtificialDiffusionStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodStrategy.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class GeometricEntity; }

      namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent an artificial diffusion strategy
/// @author Nadege Villedieu
class FluctSplit_API ArtificialDiffusionStrategy : public FluctuationSplitStrat {

  public: // typedefs

    typedef Framework::BaseMethodStrategyProvider<FluctuationSplitData,ArtificialDiffusionStrategy> PROVIDER;
    typedef Framework::VarSetMatrixTransformer MatrixTransformer;
    typedef Framework::VarSetTransformer       VectorTransformer;

  public: // methods

  /// Constructor.
  ArtificialDiffusionStrategy(const std::string& name);

  /// Destructor.
  virtual ~ArtificialDiffusionStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup the data in this object
  virtual void setup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Sets the current cell for computation of the residual splitting
  /// @param cell the pointer to the cell to treat
  void setCell(Framework::GeometricEntity* cell)
  {
      m_cell = cell;
  }

  /// Gets the current cell for computation of the residual splitting
  /// @return the pointer to the cell to treat
  Framework::GeometricEntity* getCell() const {  return m_cell;  }

  virtual void addArtificialDiff(std::vector<RealVector>& residual) = 0;

  CFuint size() const { return m_result.size(); };

  /// Gets the Class name
  static std::string getClassName() { return "ArtificialDiffusionStrategy";  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  void setIdx(const CFuint& value)  { m_idx = value; };

protected: // methods

  /// Compute the transformed states
  /// @pre assumes m_cell has been set to the current geometric entity
  /// @param states the states to be transformed to the consistent states
  /// @return the transformed consistent states
  std::vector<Framework::State*> * computeConsistentStates(std::vector<Framework::State*> *const states);

protected: // data

  /// index of the node which residual is being computed
  CFuint m_idx;
  /// vector to return with the result
  RealVector m_result;
  /// The sockets to use in this strategy for the normals
  Framework::DataSocketSink< InwardNormalsData*> socket_normals;

  /// The sockets for the flags telling if states are on the boundary
  Framework::DataSocketSink<bool> socket_isBState;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;

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

  /// the Jacobian Linearizer
  Common::SafePtr<Framework::JacobianLinearizer> _linearizer;

  /// the distribution Variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _distribVar;

  /// temporary storage of the source term
  RealVector _phiS;

  /// nodal source terms in the current cell
  std::vector<RealVector> m_nodalST;

  /// Pointer to the current cell being processed.
  Framework::GeometricEntity * m_cell;
  /// all states in current cell
  std::vector<Framework::State*> * m_states;

  /// transformed states in current cell
  std::vector<Framework::State*> * m_tStates;

  /// transformed linear states in current cell
  std::vector<Framework::State*> * m_linearStates;

  /// vector storing the flag telling if states in a cell are on the boundary
  std::vector<bool> _isBoundaryState;

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;

  /// index of the current cell
  CFuint m_cellID;

}; // class

//////////////////////////////////////////////////////////////////////////////

inline std::vector<Framework::State*> *
ArtificialDiffusionStrategy::computeConsistentStates(std::vector<Framework::State*> *const states)
{
  using namespace std;
  using namespace COOLFluiD::Framework;

  cf_assert(_updateToLinearVecTrans.isNotNull());
  m_linearStates = _updateToLinearVecTrans->transform(states);

  // computes an average state in which jacobians will be linearized
  _linearizer->linearize(*m_linearStates);

  // transformation from linear to consistent variables
  return _linearToDistMatTrans->transformFromRef(m_linearStates);
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluid

#endif // COOLFluiD_Numerics_FluctSplit_ArtificialDiffusionStrategy_hh
