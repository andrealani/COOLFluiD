#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_StdSourceTerm_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_StdSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class BlockAccumulator;
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * A base command for adding a source term
 *
 * @author Ray Vandenhoeck
 *
 *
 */
class StdSourceTerm : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit StdSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSourceTerm();
  
  /// Defines the Config Options of this class
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSockets that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * get data required for source term computation
   */
  virtual void getSourceTermData();

  /**
   * add the source term
   */
  virtual void addSourceTerm(RealVector& resUpdates) = 0;
  
  /**
   * add the src term jacobian
   */
  virtual void addSrcTermJacob();
  
  /**
   * add the src term analytical jacobian
   */
  virtual void addSrcTermJacobAna();
  
  /**
   * add the src term to the RHS
   */
  void updateRHS();
  
  virtual void getSToStateJacobian(const CFuint iState){};
  
  virtual void getSToGradJacobian(const CFuint iState){};
  
  virtual bool isGradDependent(){return false;};
  
  /// compute the gradient to state jacobian analytically
  virtual void computeGradToStateJacobianAna();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// current cell
  Framework::GeometricEntity* m_cell;

  /// cell states
  std::vector<Framework::State*>* m_cellStates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// element type index
  CFuint m_iElemType;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// solution point Jacobian determinants
  std::valarray<CFreal> m_solPntJacobDets;
  
  /// flag telling whether to add the jacobian
  bool m_addJacob;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;
  
  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;
  
  /// unperturbed updates to the residuals
  RealVector m_resUpdates;
  
  /// derivative of update to one element-residual
  RealVector m_derivResUpdates;
  
  /// number of solution points
  CFuint m_nbrSolPnts;
  
  /// perturbed solution point
  CFuint m_pertSol;
  
  /// bool telling whether the state is perturbed
  bool m_isPerturbed;
  
  /// bool telling whether analytical jacobian should be used.
  bool m_useAnaJacob;
  
  /// jacobian of the ST to the state in a solution point
  std::vector< RealVector > m_stateJacobian;

}; // class StdSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_StdSourceTerm_hh

