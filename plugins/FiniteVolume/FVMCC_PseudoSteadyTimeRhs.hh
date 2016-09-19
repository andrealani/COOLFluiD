#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhs_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_StdComputeTimeRhs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

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
class FVMCC_PseudoSteadyTimeRhs : public FVMCC_StdComputeTimeRhs {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit FVMCC_PseudoSteadyTimeRhs(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_PseudoSteadyTimeRhs();

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
  
protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;
  
  /// storage of the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;
  
  /// storage of volumes
  Framework::DataSocketSink< CFreal> socket_volumes;
  
  /// vector transformation from the update variables to the
  /// solution variables
  Common::SafePtr<Framework::VarSetTransformer> _updateToSolutionVecTrans;

  /// matrix transformation from the update variables to the
  /// solution variables starting from update variables
  Common::SafePtr<Framework::VarSetMatrixTransformer> _updateToSolutionInUpdateMatTrans;

  /// vector for the temporary flux finite difference
  RealVector _fluxDiff;

  /// temporary state
  RealVector _tempState;

  /// temporary perturbed state
  RealVector _tempPertState;

  // accumulator
  std::auto_ptr<Framework::BlockAccumulator> _acc;

  /// diagonal value in the block to be inserted in the LSS matrix
  CFreal _diagValue;

  /// flag telling if to use global DT (global time stepping)
  bool _useGlobalDT;

  /// flag telling if to use analytical transformation matrix
  bool _useAnalyticalMatrix;

}; // class FVMCC_PseudoSteadyTimeRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhs_hh
