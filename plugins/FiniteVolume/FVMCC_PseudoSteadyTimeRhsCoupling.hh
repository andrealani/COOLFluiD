#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsCoupling_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_StdComputeTimeRhsCoupling.hh"

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
class FVMCC_PseudoSteadyTimeRhsCoupling : public FVMCC_StdComputeTimeRhsCoupling {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_PseudoSteadyTimeRhsCoupling(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_PseudoSteadyTimeRhsCoupling();

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
  
protected:

  /**
   * Compute the transformation matrix numerically
   */
  virtual void computeNumericalTransMatrix(const CFuint iState,
					   const CFuint iLSS);
  
  /**
   * Compute the analytical transformation matrix
   */
  virtual void computeAnalyticalTransMatrix(const CFuint iState,
					    const CFuint iLSS);
  
protected: // data

  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;

  /// storage of the rhs
  Framework::DataSocketSink< CFreal> socket_rhs;
  
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

  // accumulator for all the LSSMatrix's
  std::vector<Framework::BlockAccumulator*> _acc;
  
  /// diagonal value in the block to be inserted in the LSS matrix
  CFreal _diagValue; 
  
  /// current number of equations
  CFuint _nbEqs;
  
  /// ID of the current starting variables
  CFuint _start;
    
}; // class FVMCC_PseudoSteadyTimeRhsCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_PseudoSteadyTimeRhsCoupling_hh
