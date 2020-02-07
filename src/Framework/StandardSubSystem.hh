// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StandardSubSystem_hh
#define COOLFluiD_Framework_StandardSubSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Stopwatch.hh"
#include "Framework/SubSystem.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Method;
    class ConvergenceMethod;
    class SpaceMethod;
    class ErrorEstimatorMethod;
    class CouplerMethod;
    class MeshAdapterMethod;
    class OutputFormatter;
    class DataProcessingMethod;
    class MeshCreator;
    class PhysicalModelImpl;
    class NumericalCommand;
    class NumericalStrategy;
    class ComputeCFL;
    class LinearSystemSolver;
    class StopConditionController;
    class SubSystemStatus;
    
//////////////////////////////////////////////////////////////////////////////

/// A StandardSubSystem is a concrete implementation of the SubSystem.
/// @author Tiago Quintino
/// @author Andrea Lani
/// @author Thomas Wuilbaut
class Framework_API StandardSubSystem : public SubSystem {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// The default constructor without arguments.
  /// @see ~StandardSubSystem()
  StandardSubSystem(const std::string& name);

  /// Destructor
  /// @see StandardSubSystem()
  virtual ~StandardSubSystem();

  /// Plugs and allocates all the sockets
  /// This complements the function on the parent class.
  /// @pre configure() as been called
  void allocateAllSockets();

  /// UnPlugs and deallocates all the sockets
  /// This complements the function on the parent class.
  /// @pre configure() as been called
  void deallocateAllSockets();

  /// Building MeshData Phase.
  void buildMeshData();
  
  /// Building Physical Model Phase.
  virtual void buildPhysicalModel();
  
  /// Setup Phase. Setup parameters of the StandardSubSystem.
  /// @see MeshData
  virtual void setup();

  /// Run (Process) Phase. All the big number crunching work goes here.
  virtual void run();

  /// Unsetup Phase.
  virtual void unsetup();

  /// Configures this Simualtion.
  /// Sets up the data for this Object.
  virtual void configure( Config::ConfigArgs& args );

  /// Adds the ActionListener's of this EventListener to the EventHandler
  void registActionListeners();

protected: // functions
  
  /// Set up MeshData that are global across partitions
  void setGlobalMeshData();
  
  /// Update the Convergence History
  void updateConvergenceFile(const CFuint nbIter);

  /// Write on the Solution on the disk
  void writeSolution(const bool forceWriting);

  /// write convergence information to stdout
  void writeConvergenceOnScreen();

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_MODIFYRESTART" Event
  /// @param eModifyRestart the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t modifyRestartAction(Common::Signal::arg_t eModifyRestart);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_AFTERGLOBALREMESHING" Event
  /// @param eAfterRemesh the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t afterRemeshingAction(Common::Signal::arg_t eAfterRemesh);

  /// Gets a vector with all the Methods this StandardSubSystem has
  /// @return vector with Method pointers.
  virtual std::vector<Framework::Method*> getMethodList();
  
  /// Set the Method's collaborators.
  virtual void setCollaborators();

  /// Set the physical model
  void configurePhysicalModel ( Config::ConfigArgs& args );

  /// Configures the coupler method tuple (some special thing need to be done here!)
  void configureCouplerMethod( Config::ConfigArgs& args, MultiMethodTuple<Framework::CouplerMethod>& couplerMethod);

  /// Sets Up the coupler method tuple (some special thing need to be done here!)
  void setCouplerMethod();

  /// Set the non root flags in the methods
  void setNonRootMethods();
  
  /// Dump the states to file
  void dumpStates();
  
  /// setup all physical models in the different namespaces
  void setupPhysicalModels();
 
  /// tell if to keep on iterating
  bool iterate(Common::SafePtr<Framework::SubSystemStatus> currSSS);
  
 protected: // data
    
  /// Duration of the simulation of this SubSystem
  Common::HourMinSec m_duration;

  /// Array of flags to be sent
  std::vector<int> m_sendFlags;
  
  /// Array of flags to be received
  std::vector<int> m_recvFlags;
  
  /// the simulation will be a is restart?
  bool m_restart;

  /// Name of Pilot file
  std::string m_pilotFileStr;

  /// PhysicalModel string
  std::string m_physicalMdlStr;

  /// configuration string for the StopCondition
  std::string m_stopConditionStr;

  /// subsystem status to be used for stop condition  
  std::string m_stopConditionSSS;

  /// Creator of the MeshData
  MultiMethodTuple<MeshCreator> m_meshCreator;

  /// MeshAdapterMethod to adapt the mesh
  MultiMethodTuple<MeshAdapterMethod> m_meshAdapterMethod;

  /// CouplerMethod to use to couple with other subsystems
  MultiMethodTuple<CouplerMethod> m_couplerMethod;

  /// ErrorEstimatorMethod to discretize the domain
  MultiMethodTuple<ErrorEstimatorMethod> m_errorEstimatorMethod;

  /// SpaceMethod to discretize the domain
  MultiMethodTuple<SpaceMethod> m_spaceMethod;

  /// ConvergenceMethod to discretize the domain
  MultiMethodTuple<ConvergenceMethod> m_convergenceMethod;

  /// LinearSystemSolver to solve linear system
  MultiMethodTuple<LinearSystemSolver> m_linearSystemSolver;

  /// DataProcessingMethod to pre process the data
  MultiMethodTuple<DataProcessingMethod> m_dataPreProcessing;

  /// DataProcessingMethod to post process the data
  MultiMethodTuple<DataProcessingMethod> m_dataPostProcessing;

  /// Output formatter
  MultiMethodTuple<OutputFormatter> m_outputFormat;
  
  /// The StopConditionController
  std::auto_ptr<StopConditionController> m_stopCondControler;

  /// Initial Physical Time Step
  CFreal m_initialTime;

  /// Initial Iteration Number
  CFuint m_initialIter;

  ///flag to force stopping the run()
  int m_forcedStop;

}; // class StandardSubSystem

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_StandardSubSystem_hh
