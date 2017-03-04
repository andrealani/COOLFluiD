// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SubSystem_hh
#define COOLFluiD_Framework_SubSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Config/ConfigObject.hh"

#include "Common/BadValueException.hh"
#include "Common/SelfRegistPtr.hh"
#include "Common/OwnedObject.hh"
#include "Common/DynamicObject.hh"

#include "Common/CFLog.hh"
#include "Environment/ConcreteProvider.hh"

#include "Framework/MultiMethodTuple.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Method;
    class NumericalCommand;
    class NumericalStrategy;
    class Namespace;
    class InteractiveParamReader;

//////////////////////////////////////////////////////////////////////////////

/// A SubSystem is configured and controlled by a Simulator.
/// This class represents a SubSystem. It is configured by Simulator through
/// feeding a startfile where the parameters are read from.
/// It should be preProcessed before processing.
/// PostProcessing is optional.
/// Memory usage function was friendly provided by Aristotelis Athanasiadis.
/// @see Simulator
/// @author Tiago Quintino
/// @author Andrea Lani
/// @author Thomas Wuilbaut
class Framework_API SubSystem :
    public Config::ConfigObject,
    public Common::DynamicObject,
    public Common::OwnedObject
{

public: // typedefs

  typedef Environment::ConcreteProvider<SubSystem,1> PROVIDER;
  typedef const std::string& ARG1;

  typedef std::vector<Common::SafePtr<NumericalCommand> >  NCommandVec_t;
  typedef std::vector<Common::SafePtr<NumericalStrategy> > NStrategyVec_t;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// The default constructor without arguments.
  SubSystem(const std::string& name);

  /// Destructor
  virtual ~SubSystem();

  /// Setup Phase. Setup parameters of the SubSystem.
  /// @pre buildMeshData() as been called
  virtual void setup() = 0;

  /// Plugs and allocates all the sockets
  /// All connections of sockets between Methods themselves,
  /// and Methods and MeshData's are established.
  /// @pre configure() as been called
  virtual void allocateAllSockets();

  /// UnPlugs and deallocates all the sockets
  /// All connections of sockets between Methods themselves,
  /// and Methods and MeshData's are broken.
  /// @pre allocateAllSockets() as been called
  virtual void deallocateAllSockets();

  /// Building MeshData Phase.
  /// @see MeshData
  /// @pre allocateAllSockets() has been called
  virtual void buildMeshData() = 0;
  
  /// Building Physical Model Phase.
  /// @see PhysicalModel
  /// @pre allocateAllSockets() has been called
  virtual void buildPhysicalModel() = 0;
  
  /// Run / Process Phase.
  /// All the big number crunching work goes here.
  /// @pre setup() as been called
  virtual void run() = 0;

  /// Unsetup Phase.
  /// @pre setup() as been called
  virtual void unsetup() = 0;

  /// Configures this Simualtion.
  /// Sets up the data for this Object.
  virtual void configure ( Config::ConfigArgs& args );

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Gets the Class name
  static std::string getClassName() {   return "SubSystem"; }

  /// Access the interactive parameter reader for the SubSystem
  Common::SafePtr<InteractiveParamReader> getInteractiveParamReader() const
  {
    cf_assert(m_int_param_reader != CFNULL);
    return m_int_param_reader;
  }

protected: // helper functions

  /// Allow (or not) to use the namespace stacks
  /// @param enable flag saying if stack usage is allowed
  void setEnableNamespaces(const bool enable);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_PLUGSOCKETS" Event
  /// @param eSocketsPlug the event which provoked this action
  /// @return an Event with a message in its body
  virtual Common::Signal::return_t allocateAllSocketsAction(Common::Signal::arg_t eSocketsAlloc);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_UNPLUGSOCKETS" Event
  /// @param eSocketsUnPlug the event which provoked this action
  /// @return an Event with a message in its body
  virtual Common::Signal::return_t deallocateAllSocketsAction(Common::Signal::arg_t eSocketsDealloc);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_BUILDMESHDATA" Event
  /// @param eSetup the event which provoked this action
  /// @return an Event with a message in its body
  virtual Common::Signal::return_t buildMeshDataAction(Common::Signal::arg_t eBuild);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_BUILDPHYSICALMODEL" Event
  /// @param eSetup the event which provoked this action
  /// @return an Event with a message in its body
  virtual Common::Signal::return_t buildPhysicalModelAction(Common::Signal::arg_t eBuild);
  
  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_SETUP" Event
  /// @param eSetup the event which provoked this action
  /// @return an Event with a message in its body
  virtual Common::Signal::return_t setupAction(Common::Signal::arg_t eSetup);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_RUN" Event
  /// @param eRun the event which provoked this action
  /// @return an Event with a message in its body
  virtual Common::Signal::return_t runAction(Common::Signal::arg_t eRun);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_UNSETUP" Event
  /// @param eUnsetup the event which provoked this action
  /// @return an Event with a message in its body
  virtual Common::Signal::return_t unsetupAction(Common::Signal::arg_t eUnsetup);

  /// Gets all the commands in all the Methods
  std::pair < NCommandVec_t, NStrategyVec_t >  getAllCommandsAndStrategies(); // const;

  /// Setups all the Namespace's each with a MeshData,
  /// PhysicalModel and a SubSystemStatus
  /// If none have been configuredm then create a "Default" namespace with a
  /// Default MeshData / PhysicalModel and SubSystemStatus
  void configureNamespaces(Config::ConfigArgs& args);

  /// Creates the Singletons on their NamespaceStack's
  /// This currently applies to MeshData, PhysicalModel and SubSystemStatus
  void configureSingletons( Config::ConfigArgs& args );

  /// Alloocate and configure the Singletons if they do not exist yet
  void configureNamespaceSingletons( Config::ConfigArgs& args, Common::SafePtr<Namespace> nsp);

  /// Set all NumericalCommand's to prepare the SubSystem
  void setCommands();

  /// Puts all the commands from the Method
  void putCommandsFromMethod(Method* mtd, NCommandVec_t& coms);// const;

  /// Puts all the strategies from the Method
  void putStrategiesFromMethod(Method* mtd, NStrategyVec_t& strat);// const;

  /// Template function for creating and configuring multi Method's
  template <class METHOD>
  void configureMultiMethod(Config::ConfigArgs& args, MultiMethodTuple<METHOD>& method);

  /// Template function for creating and configuring Method's
  template <class BASEMETHOD>
  void configureMethod(Config::ConfigArgs& args,
                       Common::SelfRegistPtr<BASEMETHOD>& mtd,
                       const std::string& type,
                       const std::string& name = "");

  /// Set the collaborators in the corresponding to the methods.
  template <class COLLAB>
  void setMatchingCollaborators(const std::vector<QualifiedName>& names,
                                MultiMethodTuple<COLLAB>& collab,
                                MultiMethodHandle<COLLAB>& cList);

  /// Gets a vector with all the Methods this SubSystem has
  /// @return vector with Method pointers.
  virtual std::vector<Framework::Method*> getMethodList() = 0;
  
  /// Allocates all sockets in all methods of this SubSystem
  void allocateSockets();

  /// Deallocates all sockets in all methods of this SubSystem
  void deallocateSockets();

  /// Matches all sockets in all methods of this SubSystem
  void matchAllSockets();

  /// Check that all essential DataSocketSink's in the Method's
  /// of this SubSystem have been satisfied.
  void checkAllSockets();

  /// Sets the Parent Namespace in the Methods Sockets
  void setParentNamespaceInMethodSockets();
  
  /// fill in the array holding group ranks 
  void fillGroupRanks(const std::string rankString, 
		      std::vector<int>& granks);
  
  /// Set the Method's collaborators.
  virtual void setCollaborators() = 0;

  /// Set the Method's collaborators.
  template <class METHOD, class COLLAB>
  void setCollaborators(MultiMethodTuple<METHOD>& method,
                        MultiMethodTuple<COLLAB>& collaborator);
  
  /// @return the counter of the active ranks involved in this subsystem
  CFuint getNbActiveRanks() const {return m_ranksCounter;}
  
 private: // member data

  /// counter for the number of ranks to be activated in this SubSystem
  CFuint m_ranksCounter;
  
  /// names of the namespaces of this SubSystem
  std::vector<std::string> m_namespaces;
  
  /// ranks associated to each namespace in the form
  /// START0:END0 START1:END1 etc.
  std::vector<std::string> m_ranks;
  
  /// indicates to create null methods
  bool m_has_null_methods;

  /// number of methods allocated so far.
  /// This is used to create unique names for the null methods.
  CFuint m_nbmethods;

  /// Interactive parameter reader
  InteractiveParamReader * m_int_param_reader;

}; // class SubSystem

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(SubSystem)

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystem.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SubSystem_hh
