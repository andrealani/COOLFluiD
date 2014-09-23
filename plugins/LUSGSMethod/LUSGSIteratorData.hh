#ifndef COOLFluiD_Numerics_LUSGSMethod_LUSGSIteratorData_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUSGSIteratorData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeNorm.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/NumericalJacobian.hh"
#include "Framework/SubSystemStatus.hh"

#include "LUSGSMethod/ComputeNormLUSGS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * LUSGSMethodCom 's that compose the LUSGSMethod.
   *
   * @see LUSGSMethodCom
   *
   * @author Kris Van den Abeele
   * @author Matteo Parsani
   */
class LUSGSIteratorData : public Framework::ConvergenceMethodData {

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  LUSGSIteratorData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~LUSGSIteratorData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "LUSGSIterator";
  }

  Common::SafePtr<ComputeNormLUSGS> getLUSGSNormComputer() const
  {
    return m_computeNormLUSGS;
  }

  /**
   * Get the numerical jacobian calculator
   */
  Common::SafePtr<Framework::NumericalJacobian> getNumericalJacobian() const
  {
    cf_assert(m_numJacob.get() != CFNULL);
    return m_numJacob.get();
  }

  /**
   * Gets the number of time-steps to perform before recomputing the block Jacobiam matrices
   * @param iSys is the ID of the current subsystem of equations to solve
   */
  CFuint getJacobFreezFreq(CFuint iSys) const
  {
    cf_assert(iSys < m_jacobFreezFreq.size());
    return m_jacobFreezFreq[iSys];
  }

  /**
   * Gets the maximum number of sweeps to be performed
   * @param iSys is the ID of the current subsystem of equations to solve
   */
  CFuint getNbMaxSweeps(CFuint iSys) const
  {
    cf_assert(iSys < m_maxSweeps.size());
    return m_maxSweeps[iSys];
  }

  /**
   * Gets the maximum norm to be achieved
   */
  CFreal getMaxNorm() const
  {
    return m_maxNorm;
  }

  /**
   * Checks if convergence history should be printed.
   */
  bool isPrintHistory() const
  {
    return m_printHistory;
  }

  /**
   * Gets the flag that indicates we are at the last iteration
   */
  bool isAchieved() const
  {
    return m_achieved;
  }

  /**
   * Sets the flag that indicates we are at the last iteration
   */
  void setAchieved(bool achieved)
  {
    m_achieved = achieved;
  }

  /**
   * Gets the flag that indicates we are at the forward sweep.
   * @returns m_forwardSweep
   */
  bool isForwardSweep()
  {
    return m_forwardSweep;
  }

  /**
   * Sets m_forwardSweep
   */
  void setForwardSweep(const bool forwardSweep)
  {
    m_forwardSweep = forwardSweep;
  }

  /**
   * Gets the flag that indicates we are at the last states set.
   * @returns m_stopSweep
   */
  bool stopSweep()
  {
    return m_stopSweep;
  }

  /**
   * @returns m_stopStatesLoop
   */
  bool stopStatesLoop()
  {
    return m_stopStatesLoop;
  }

  /**
   * @returns m_stopEqsLoop
   */
  bool stopEqsLoop()
  {
    return m_stopEqsLoop;
  }

  /**
   * @returns m_beforePertResComputation
   */
  bool beforePertResComputation()
  {
    return m_beforePertResComputation;
  }

  /**
   * Sets m_stopSweep
   */
  void setStopSweep(const bool stopSweep)
  {
    m_stopSweep = stopSweep;
  }

  /**
   * Sets m_stopStatesLoop
   */
  void setStopStatesLoop(const bool stopStatesLoop)
  {
    m_stopStatesLoop = stopStatesLoop;
  }

  /**
   * Sets m_stopEqsLoop
   */
  void setStopEqsLoop(const bool stopEqsLoop)
  {
    m_stopEqsLoop = stopEqsLoop;
  }

  /**
   * Sets m_beforePertResComputation
   */
  void setBeforePertResComputation(const bool beforePertResComputation)
  {
    m_beforePertResComputation = beforePertResComputation;
  }

  /**
   * Gets the flag that indicates the number of states sets minus one.
   * @return m_nbrStatesSetsMinusOne
  */
  CFuint getNbrStatesSets()
  {
    return m_nbrStatesSets;
  }
  /**
   * Sets the flag that indicates the number of states sets minus one.
   * @return m_nbrStatesSetsMinusOne
   */
  void setNbrStatesSets(const CFuint nbrStatesSets)
  {
    m_nbrStatesSets = nbrStatesSets;
  }

  /**
   * @return a pointer to m_resAux
   */
  Common::SafePtr< RealVector > getResAux()
  {
    return &m_resAux;
  }

  /**
   * @return m_withPivot
   */
  bool isWithPivot()
  {
    return m_withPivot;
  }

  /**
   * @return m_withPivot
   */
  void setIsWithPivot(const bool withPivot)
  {
    m_withPivot = withPivot;
  }

private: // data

  /// Functor that computes the requested norm specific for LUSGSMethod
  Common::SafePtr<ComputeNormLUSGS>  m_computeNormLUSGS;

  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> m_numJacob;

  /// flag to indicate that convergence has been achieved
  bool m_achieved;

  /// number of time-steps to perform before recomputing the block Jacobiam matrices
  std::vector<CFuint> m_jacobFreezFreq;

  /// number of sweeps per lusgs iterations
  std::vector<CFuint> m_maxSweeps;

  /// L2 norm of dU to reach per global time step
  CFreal m_maxNorm;

  /// flag to indicate printing of history in the newton iteration
  bool m_printHistory;

  /// boolean telling whether it is a forward or a backward sweep
  bool m_forwardSweep;

  /// boolean telling whether last statesSet is reached
  bool m_stopSweep;

  /// boolean telling whether last state in current statesSet is reached
  bool m_stopStatesLoop;

  /// boolean telling whether last variable in current state
  bool m_stopEqsLoop;

  /// boolean telling whether we are before the computation of the perturbed residual
  bool m_beforePertResComputation;

  /// number of states sets
  CFuint m_nbrStatesSets;

  /// Residual auxiliary variable
  RealVector m_resAux;

  /// boolean telling whether pivotation is used
  bool m_withPivot;

}; // end of class LUSGSIteratorData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for LUSGSMethod
typedef Framework::MethodCommand<LUSGSIteratorData> LUSGSIteratorCom;

/// Definition of a command provider for LUSGSMethod
typedef Framework::MethodCommand<LUSGSIteratorData>::PROVIDER LUSGSIteratorComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_LUSGSIteratorData_hh
