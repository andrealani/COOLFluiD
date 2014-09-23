#ifndef COOLFluiD_Numerics_FluctSplit_ComputeRHS_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeRHS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

      class FluctuationSplitStrategy;

//////////////////////////////////////////////////////////////////////////////

/// This class represent a command that computes the RHS
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API ComputeRHS : public FluctuationSplitCom {
public:

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer VectorTransformer;

  /// Constructor.
  explicit ComputeRHS(const std::string& name);

  /// Destructor.
  ~ComputeRHS();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// UnSet up private data and data of the aggregated classes
  /// in this command after processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // methods

  /// Execute the command on the current TRS
  virtual void executeOnTrs();

  /// Cleans the rhs setting it to zero
  virtual void cleanRHS();

  /// Transform the residual
  void transformResidual();

protected: // data

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

  /// Matrix Transformer from Solution to Update in Update Variables
  Common::SafePtr<MatrixTransformer> _solutionToUpdateMatTrans;

  /// Transformer from Update to Linear Variables
  Common::SafePtr<VectorTransformer> _updateToLinearVecTrans;

  /// the distribution Variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _distribVar;

  /// the update Variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _updateVar;

  /// the linearization Variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _linearVar;

  /// diffusive term computer
  Common::SafePtr<ComputeDiffusiveTerm> _diffTermComputer;

  /// socket for storage of the rhs
  Framework::DataSocketSink< CFreal> socket_rhs;

  /// socket for storage of the update coefficients
  Framework::DataSocketSink< CFreal> socket_updateCoeff;

  /// socket for the inward normals storage
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// socket for the flags telling if states are on the boundary
  Framework::DataSocketSink<bool> socket_isBState;

  /// storage of the States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for the flags telling if states have been updated
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// temporary storage of the cell residual
  std::vector<RealVector> _residual;

  /// temporary storage of the diffusive cell residual
  std::vector<RealVector> _diffResidual;
  /// temporary storage of the diffusive cell residual
  std::vector<RealVector> _artdiffResidual;
  /// fluctuation split strategy object
  Common::SafePtr<FluctuationSplitStrategy> _fsStrategy;
  /// fluctuation split strategy object
  Common::SafePtr<ArtificialDiffusionStrategy> _adStrategy;

  /// flag telling if a diffusive term has to be computed
  bool _hasDiffusiveTerm;

  /// flag telling if a diffusive term has to be computed
  bool _hasArtDiffusiveTerm;

  /// flag forcing to freeze the diffusive coefficients
  /// while doing numerical perturbation of the jacobians
  bool _freezeDiffCoeff;

}; // class ComputeRHS

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeRHS_hh

