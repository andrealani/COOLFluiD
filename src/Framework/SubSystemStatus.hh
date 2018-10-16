// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SubSystemStatus_hh
#define COOLFluiD_Framework_SubSystemStatus_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Stopwatch.hh"

#include "Framework/ComputeDT.hh"
#include "Framework/ConvergenceStatus.hh"
#include "Framework/NamespaceStack.hh"
#include "Framework/VarRegistry.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Common { 
    class FactoryRegistry;
  }
  
  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

/// This class represents a singleton object where global infos (number
/// of iterations and residual) about the status of the SubSystem
/// are held
/// @author Andrea Lani
/// @author Dries Kimpe
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API SubSystemStatus : public Common::NonCopyable<SubSystemStatus>,
				      public Config::ConfigObject {

friend class SubSystemStatusStack;

public: // methods

  /// Set the factory registry
  void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr);
  
  /// Get the factory registry
  Common::SafePtr<Common::FactoryRegistry> getFactoryRegistry();
  
  /// Gets the variable registry
 Common::SafePtr<VarRegistry> getVarRegistry() {return m_var_registry;}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Get some convergence data bundled into an object
  ConvergenceStatus getConvergenceStatus();

//////////////////////////////////////////////////////////////////////////////
//
//   Residual related functions
//
//////////////////////////////////////////////////////////////////////////////

  /// Sets the flag asking the simulation to be stopped
  void setStopSimulation(const bool stopSim) {m_stopSim = stopSim;}
  
  /// @return the flag asking the simulation to be stopped
  bool getStopSimulation() const {return m_stopSim;}
    
  /// Sets the residual
  void setResidual(const RealVector& residual);
  
  /// Sets the residual
  void setMonitoredVar(const CFuint monitored_var)
  {
    m_monitored_var = monitored_var;
  }
  
  /// Set the global residual to use in the convergence method
  void setGlobalRes(const bool global_res)
  {
    m_global_res = global_res;
  }

  /// Get the value of the residual of the monitored variable
  /// @see setMonitoredVar
  /// @post if m_residual vector has no size then no ConvergenceMethod
  ///       has set the residual in the SubSystemStatus, therefore return 0
  /// @return the residual of the monitored variable
  CFreal getResidual() const;


  /// Get the value of the residual
  /// @return the residuals of all the cpomputed variables
  RealVector getAllResiduals() const { return m_residual; }


  /// Reset to zero the residual
  void resetResidual() { m_residual = 0.0; }

//////////////////////////////////////////////////////////////////////////////
//
//   Iteration related functions
//
//////////////////////////////////////////////////////////////////////////////


  /// Get the number of iterations
  CFuint getNbIter() const { return m_iter; }

  /// Reset to zero the iteration
  void setNbIter(const CFuint p) { m_iter = p; }

  /// Update the number of iterations by incrementing it
  void updateNbIter();

//////////////////////////////////////////////////////////////////////////////
//
//   Iteration/Timestep InnerStep related functions
//
//////////////////////////////////////////////////////////////////////////////

  /// Get the info if it is setup phase
  bool isSetup() const { return m_isSetup; }
  
  /// Get the info if it is setup phase
  void setSetup(bool p) { m_isSetup = p; }
  
  /// Get the info if first step
  bool isFirstStep() const { return m_firstStep; }

  /// Get the info if last step
  bool isLastStep() const { return m_lastStep; }

  /// Set the info if first step
  void setFirstStep(bool p) { m_firstStep = p; }

  /// Set the info if last step
  void setLastStep(bool p) { m_lastStep = p; }

//////////////////////////////////////////////////////////////////////////////
//
//   Moving Mesh related functions
//
//////////////////////////////////////////////////////////////////////////////

 
  
  /// Get the info if moving mesh
  bool isMovingMesh() const { return m_movingMesh; }

  /// Set the info if moving mesh
  void setMovingMesh(bool p)  { m_movingMesh = p; }

  /// Set the info to append to outputfiles
  void setAppendToFile(const bool iter, const bool time)
  {
    m_appendIter = iter;
    m_appendTime = time;
  }

  /// Needs to append the iteration number??
  bool isAppendIter() { return m_appendIter; }

//////////////////////////////////////////////////////////////////////////////
//
//   Time and Time Step related functions
//
//////////////////////////////////////////////////////////////////////////////


  /// Needs to append the subsystem Time??
  bool isAppendTime() { return m_appendTime; }

  /// Reset to zero the current time
  void resetCurrentTime() { m_currentTime = 0.; }
  
  /// Update the current Time using the current Time Step
  void updateCurrentTime();
  
  /// Get the current (adimensional) Time
  CFreal getCurrentTime() const { return m_currentTime; }

  /// Get the current dimensional Time
  CFreal getCurrentTimeDim() const;

  /// Set the current (adimensional) Time
  void setCurrentTime(CFreal const p)  { m_currentTime = p;}

  /// Set the current dimensional Time
  void setCurrentTimeDim(CFreal const p);

  /// Get the maximum expected (adimensional) time of simulation
  CFreal getMaxTime() const { return m_max_time; }

  /// Get the maximum expected dimensional time of simulation
  CFreal getMaxTimeDim() const;

  /// Set the maximum expected (adimensional) time of simulation
  void setMaxTime(CFreal time) { m_max_time = time; }

  /// Set the maximum expected dimensional time of simulation
  void setMaxTimeDim(CFreal time);

  /// Get the maximum (adimensional) TimeStep (stability)
  CFreal getMaxDT() const { return m_maxDT; }

  /// Get the maximum dimensional TimeStep (stability)
  CFreal getMaxDTDim() const;

  /// Get the maximum (adimensional) TimeStep (stability)
  void setMaxDT(CFreal p) { m_maxDT = p; }

  /// Get the maximum adimensional TimeStep (stability)
  void setMaxDTDim(CFreal p);

  /// Get the (adimensional) DT
  CFreal getDT() const  { return m_timeStep;}

  /// Get the (adimensional) DT
  CFreal getDTDim() const;

  /// Get the previous (adimensional) DT
  CFreal getPreviousDT() const
  {
    return (m_previousTimeStep != 0.) ? m_previousTimeStep : m_timeStep;
  }

  /// Get the previous previous (adimensional) DT
  CFreal getPrevPrevDT() const
  {
    return (m_prevprevTimeStep != 0.) ? m_prevprevTimeStep : getPreviousDT();
  }

  /// Set the (adimensional) DT
  void setDT(const CFreal DT)
  {
    m_prevprevTimeStep = m_previousTimeStep;
    m_previousTimeStep = m_timeStep;
    m_timeStep = DT;
  }

  /// Set the dimensional DT
  void setDTDim(const CFreal DT);

  /// Update the value of DT
  void updateTimeStep()  { (*m_computeDT)(); }

  /// Make the time data non-dimensional (to be called at the beginning of the computation
  void adimensionalizeTimeData();

//////////////////////////////////////////////////////////////////////////////
//
//   Multiple Time Step layers
//
//////////////////////////////////////////////////////////////////////////////


  /// Get the number of layers
  CFreal getTimeStepLayers() const { return m_timeStepLayers; }

  /// Set the inner DT
  void setInnerDTSize(const CFuint i)
  {
    m_innerDT.resize(i);
  }

  /// Set the inner DT
  void setInnerDT(const CFuint i, const CFreal DT)
  {
    cf_assert(i < m_innerDT.size());
    m_innerDT[i] = DT;
  }

  /// Set the inner DT Ratio
  void setInnerDTRatio(const CFuint i, const CFreal DTRatio)
  {
    cf_assert(i < m_innerDT.size());
    m_innerDTRatio[i] = DTRatio;
  }

  /// Get the inner DT
  CFreal getInnerDT(const CFuint i) const
  {
    cf_assert(i < m_innerDT.size());
    return m_innerDT[i];
  }

  /// Get the inner DT Ratio
  CFreal getInnerDTRatio(const CFuint i) const
  {
    cf_assert(i < m_innerDTRatio.size());
    return m_innerDTRatio[i];
  }

//////////////////////////////////////////////////////////////////////////////
//
//   SubIterations
//
//////////////////////////////////////////////////////////////////////////////


  /// Set the flag telling if we are doing subiterations
  void setSubIterationFlag(bool flag)
  {
    m_subIterationFlag = flag;
  }

  /// Set the flag telling if we are doing subiterations
  bool doingSubIterations()
  {
    return m_subIterationFlag;
  }

  /// Get the flag telling if during the first subiter of a subiteration
  /// or if we are not using subiterations
  bool isSubIterationFirstStep()
  {
    if (m_subIterationFlag && m_subIter == 0) return true;
    if (!m_subIterationFlag) return true;
    return false;
  }

  /// Get the flag telling if we are in the last subiter of a subiteration
  /// or if we are not using subiterations
  void setIsSubIterationLastStep(bool flag)
  {
    m_lastSubIter = flag;
  }

  /// Get the flag telling if we are in the last subiter of a subiteration
  /// or if we are not using subiterations
  bool isSubIterationLastStep()
  {
    if(m_subIterationFlag) return m_lastSubIter;

    return true;
  }

  /// Set the subIteration counter
  void setSubIter(const CFuint p)
  {
    m_subIter = p;
  }

  /// Get the subIteration counter
  CFuint getSubIter()
  {
    return m_subIter;
  }

  /// Update the number of sub-iterations by incrementing it
  void updateSubIter()
  {
    ++m_subIter;
  }

//////////////////////////////////////////////////////////////////////////////
//
//   Timing related functions
//
//////////////////////////////////////////////////////////////////////////////

  /// Start stopwatch
  void startWatch() { m_stopwatch.start(); }

  /// Stop stopwatch
  void stopWatch()  { m_stopwatch.stop();  }

  /// Read stopwatch
  CFreal readWatch() { return m_stopwatch.read(); }

  /// Read stopwatch using HMS format
  Common::HourMinSec readWatchHMS() { return m_stopwatch.readTimeHMS(); }

//////////////////////////////////////////////////////////////////////////////
//
//   General functions
//
//////////////////////////////////////////////////////////////////////////////


  /// Configure
  virtual void configure ( Config::ConfigArgs& args );

  /// Get the sub system name
  std::string getSubSystemName() const
  {
    return m_subSystemName;
  }

  /// Set the sub system name
  void setSubSystemName(std::string const name)
  {
    m_subSystemName = name;
  }

  /// Default destructor
  ~SubSystemStatus();

private: // methods

  /// Constructor
  SubSystemStatus(const std::string& name);

private: // member data
  
  /// factory registry to allow polymorphic creation of objects
  Common::SafePtr<Common::FactoryRegistry> m_fr;
  
  /// registry for dynamic created variables
  VarRegistry * m_var_registry;

  /// Cronometer to time the simulation
  Common::Stopwatch<Common::CPUTime> m_stopwatch;

  /// the number of iterations
  CFuint m_iter;

  /// the number of subIterations
  CFuint m_subIter;

  /// the current time (unsteady computations)
  CFreal m_currentTime;

  /// the Simulation residual number stored
  RealVector m_residual;

  /// Maximum Allowable Time Step
  CFreal m_maxDT;
  
  /// the index of the variable for which the residual is monitored
  CFuint m_monitored_var;

  /// Flag to set the global residual
  bool m_global_res;

  ///flag: is it the first step (not time step but for example newton step)
  bool m_firstStep;

  ///flag: is it the last step (not time step but for example newton step)
  bool m_lastStep;

  ///flag: tells if it is a setup phase
  bool m_isSetup;
  
  ///flag: is this subsystem using a moving mesh
  bool m_movingMesh;
  
  ///flag: should the iteration number be append to the filename
  bool m_appendIter;

  ///flag: should the current time be append to the filename
  bool m_appendTime;

  ///Do we do subiterations
  bool m_subIterationFlag;
  
  ///Is it the last step of the subiteration process
  bool m_lastSubIter;
  
  ///flag: tells the simulation to stop
  bool m_stopSim;
  
  /// Name of the subsystem
  std::string m_subSystemName;

  /// vector of the inner DeltaT
  /// used for computation
  RealVector m_innerDT;

  /// vector of the inner DeltaT
  /// used for computation
  RealVector m_innerDTRatio;

  /// vector of the inner DeltaT
  /// used at configuration time
  std::vector<CFreal> m_innerDTConf;

  /// the value of the previous DT
  CFreal m_previousTimeStep;
  
  /// the value of the previous previous DT
  CFreal m_prevprevTimeStep;
    
  /// maximum time expected time for simulation
  CFreal m_max_time;
  
  /// the TimeStep value
  CFreal m_timeStep;
  
  /// Number of DT to be stored (for multi-layer spacetime,for example)
  CFuint m_timeStepLayers;
  
  /// string for configuration of the DT term computer
  std::string m_computeDTStr;

  /// DT calculator
  Common::SelfRegistPtr<ComputeDT> m_computeDT;
  
}; // end of class SubSystemStatus

//////////////////////////////////////////////////////////////////////////////

class Framework_API SubSystemStatusStack : public NamespaceStack<SubSystemStatus> {
public:


  /// Returns the instance of the Active SubSystemStatus
  /// which is the one on top of the stack
  /// @return SafePtr to the active SubSystemStatus

  static Common::SafePtr<SubSystemStatus> getActive();


  /// Returns the instance of this meshDataStack
  /// This is the access point to the Singleton
  /// @return the instance of the singleton
  
  static SubSystemStatusStack& getInstance();
  
  /// Set the current subsystem name
  static void setCurrentName(const std::string& name) {getInstance().m_ssName = name;}
  
  /// Returns the current subsystem name
  static std::string getCurrentName() {return getInstance().m_ssName;}
  
protected: // helper functions from NamespaceStack
  
  /// Gets the name of the SubSystemStatus from the Namespace
  /// @param nsp the Namespace from where to get te object name

  std::string getObjectName(const Common::SafePtr<Namespace>& nsp) const;


  /// Creates a SubSystemStatus with the supplied name
  /// @param name of the SubSystemStatus

  SubSystemStatus * createObject(const std::string& name);

private:

  /// name of the current subsystem
  std::string m_ssName;
  
}; // end of class SubSystemStatusStack;

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SubSystemStatus_hh
