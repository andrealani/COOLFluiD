#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsBlockDiag_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsBlockDiag_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_StdComputeTimeRhs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {class BlockAccumulatorBase;}

  namespace MathTools {class MatrixInverter;}
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the pseudo steady RHS
 * using standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_PseudoSteadyTimeRhsBlockDiag : public FVMCC_StdComputeTimeRhs {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit FVMCC_PseudoSteadyTimeRhsBlockDiag(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_PseudoSteadyTimeRhsBlockDiag();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Execute Processing actions
   */
  virtual void execute();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private:

  /**
   * Compute the transformation matrix numerically
   */
  virtual void computeNumericalTransMatrix(const CFuint iState);
  
  /**
   * Compute the analytical transformation matrix
   */
  virtual void computeAnalyticalTransMatrix(const CFuint iState);

  /**
   * Solve the linear system
   */
  virtual void solveSys();
  
protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;
  
  /// storage of the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;
  
  /// storage of volumes
  Framework::DataSocketSink< CFreal> socket_volumes;

  /// storage of diagonal block matrix
  Framework::DataSocketSink<CFreal> socket_diagMatrix;
  
  /// vector transformation from the update variables to the
  /// solution variables
  Common::SafePtr<Framework::VarSetTransformer> _updateToSolutionVecTrans;

  /// matrix transformation from the update variables to the
  /// solution variables starting from update variables
  Common::SafePtr<Framework::VarSetMatrixTransformer> _updateToSolutionInUpdateMatTrans;
   
  /// temporary data for holding the matrix inverter
  std::auto_ptr<MathTools::MatrixInverter>  _inverter;
  
  /// vector for the temporary flux finite difference
  RealVector _fluxDiff;

  /// temporary state
  RealVector _tempState;

  /// temporary perturbed state
  RealVector _tempPertState;

  /// storage for a local rhs array
  std::vector<CFreal> _tempRHS;
  
  /// block accumulator
  std::auto_ptr<Framework::BlockAccumulatorBase> _acc;
  
  /// diagonal value in the block to be inserted in the LSS matrix
  CFreal _diagValue;

  /// flag telling if to use global DT (global time stepping)
  bool _useGlobalDT;

  /// flag telling if to use analytical transformation matrix
  bool _useAnalyticalMatrix;

}; // class FVMCC_PseudoSteadyTimeRhsBlockDiag

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsBlockDiag_hh
