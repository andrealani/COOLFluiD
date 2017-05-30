// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Method_hh
#define COOLFluiD_Framework_Method_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"

#include "Common/OwnedObject.hh"
#include "Common/DynamicObject.hh"
#include "Common/SetupObject.hh"

#include "Common/DynamicFunctionCaller.hh"
#include "Common/NullableObject.hh"

#include "Common/CFLog.hh"
#include "Environment/Factory.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/QualifiedName.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class MethodData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Method.
/// A Method is composed of several NumericalCommand's
/// The NumericalCommand's can be configured and interchanged.
/// @see ConfigObject
/// @see NumericalCommands
/// @author Tiago Quintino

class Framework_API Method :

    public Common::NonCopyable<Method>,
    public Common::OwnedObject,
    public Config::ConfigObject,
    public Common::SetupObject,
    public Common::DynamicObject,
    public Common::NullableObject,
    public Framework::NamespaceMember {

public: // functions

  /// Defines the user configurable options of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit Method(const std::string& name);

  /// Virtual destructor
  /// Calls unsetup() if it is still setup
  virtual ~Method();

  /// Configure this object with user defined parameters
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// @return the qualified name of the Method
  QualifiedName QName() const;

  /// Sets up the data for this method to be applied
  /// @post pushs and pops the Namespace to which this Method belongs
  void setMethod();

  /// Clears the data for this method.
  /// @post pushs and pops the Namespace to which this Method belongs
  void unsetMethod();

  /// Gets a vector with all the commands this method will use.
  /// @return vector with NumericalCommand pointers.
  virtual std::vector< Common::SafePtr<NumericalCommand> > getCommandList() const;

  /// Execute a given sequence of commands identified by name
  void executeCommands(const std::vector<std::string>& comNames);

  /// Gets a vector with all the strategies this method will use.
  /// @return vector with NumricalStrategy pointers.
  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList() const;

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr<Framework::MethodData> getMethodData() const = 0;

  /// Sets the command groups on the commands that belong to them
  void setCommandGroups();

  /// Gets the all the command groups on this Method
  std::vector<Common::SafePtr<CommandGroup> > getCommandGroups();

  /// Gets the number of groups
  CFuint getNbGroups() { return m_groupNames.size(); }

  /// Sets the Namespace on all sockets of this Method.
  /// By default all sockets belong to the SocketNamespace
  /// of the Method to which they belong
  void setSocketNamespaces();

  /// Allocates all sockets in all Commands of this Method
  void allocateMethodSockets();

  /// Deallocates all sockets in all Commands of this Method
  void deallocateMethodSockets();

  /// Matches all sockets in all Commands of this Method
  void matchMethodSockets(std::vector<Method*> methodList);

  /// Get the needed sockets for this Method
  std::vector<Common::SafePtr<BaseDataSocketSink> >  getNeededSockets();

  /// Get the provided sockets for this Method
  std::vector<Common::SafePtr<BaseDataSocketSource> > getProvidedSockets();

  /// Get All DataSockets on this Method, both Sinks and Sources
  std::vector<Common::SafePtr<DataSocket> > getAllSockets();

  /// Check that all essential DataSocketSink's in this Method have been satisfied.
  void checkAllSockets();

  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with no parameters.
  virtual void run_function(const std::string & func) = 0;

  /// Get flag to know if method is to be called by a plugin and not by the subsystem
  bool isNonRootMethod() { return m_isNonRootMethod; }

  /// Set flag to know if method is to be called by a plugin and not by the subsystem
  void setIsNonRootMethod() { m_isNonRootMethod = true; }

  /// Eventually set the collaborating methods to be non-root
  virtual void setNonRootMethods() {  /* by default nothing is done here */ }

protected: // interface

  /// Sets up the data for this method to be applied
  /// This is the abstract function that the concrete methods must implement.
  virtual void setMethodImpl() = 0;

  /// Clears the data for this method.
  /// This is the abstract function that the concrete methods must implement.
  virtual void unsetMethodImpl() = 0;

protected: // helper functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Executed by the the "CF_ON_MESH_UPDATE" Event
  /// @param generic argument for signal
  /// @return generic return from signal execution
  Common::Signal::return_t
  resetup ( Common::Signal::arg_t input );

  /// Switch to the Namespace of this Method
  void pushNamespace();

  /// Switch back from the Namespace of this Method
  void popNamespace();

  /// Configures the Command Groups in this Method
  void configureCommandGroups ( Config::ConfigArgs& args );

  /// Setups all the Commands and all the Strategies
  void setupCommandsAndStrategies();

  /// Unsetups all the Commands and all the Strategies
  void unsetupCommandsAndStrategies();

  /// Template function for creating and configuring MethodCommand's passing with
  /// stronger checking of template parameters.
  /// @param com the pointer which should be holding the newly created and configured command
  /// @param type the string with the name of the type of command to create
  /// @param name the name to give to the newly created command
  /// @param data the MethodData that the command will share with other commands and strategies of the same method
  template < typename DATA, typename PROVIDER >
  void configureCommand(Config::ConfigArgs& args,
                        Common::SelfRegistPtr< MethodCommand<DATA> >& com,
                        const std::string& type,
                        const std::string& name,
                        const Common::SharedPtr<DATA>& data)
  {
    typedef MethodCommand<DATA> BASECOMMAND;
    configureCommand<BASECOMMAND,DATA,PROVIDER>(args,com,type,name,data);
  }

  /// Template function for creating and configuring MethodCommand's without giving any name.
  /// Name will be equal to the type.
  /// @param com the pointer which should be holding the newly created and configured command
  /// @param type the string with the name of the type of command to create
  /// @param data the MethodData that the command will share with other commands and strategies of the same method
  template < typename DATA, typename PROVIDER >
  void configureCommand(Config::ConfigArgs& args,
                        Common::SelfRegistPtr< MethodCommand<DATA> >& com,
                        const std::string& type,
                        const Common::SharedPtr<DATA>& data)
  {
    typedef MethodCommand<DATA> BASECOMMAND;
    configureCommand<BASECOMMAND,DATA,PROVIDER>(args,com,type,type,data);
  }

  /// Template function for creating and configuring MethodCommand's with reduced
  /// checking of template parameters because BASECOMMAND can be more specific than a MethodCommand.
  /// @param com the pointer which should be holding the newly created and configured command
  /// @param type the string with the name of the type of command to create
  /// @param name the name to give to the newly created command
  /// @param data the MethodData that the command will share with other commands and strategies of the same method
  template < typename BASECOMMAND, typename DATA, typename PROVIDER>
  void configureCommand(Config::ConfigArgs& args,
                        Common::SelfRegistPtr<BASECOMMAND>& com,
                        const std::string& type,
                        const std::string& name,
                        const Common::SharedPtr<DATA>& data)
  {
    CFLogDebugMed( "Method::configureCommand() => Type: " << type << " Name: " << name << "\n");
    Common::SafePtr<PROVIDER> prov;
    try
    {
      prov = FACTORY_T_GET_PROVIDER(this->getFactoryRegistry(), BASECOMMAND, type);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      throw;
    }
    cf_assert(prov.isNotNull());
    com = prov->create(name,data);
    cf_assert(com.isNotNull());
    com->setFactoryRegistry(this->getFactoryRegistry());
    m_commands.push_back(com.getPtr());
    configureNested ( com.getPtr(), args );
  }

protected: // data

  /// storage of the commands in this method
  std::vector < Common::SafePtr<NumericalCommand> > m_commands;

  /// Names of the CommandGroup's created
  std::vector<std::string> m_groupNames;

  /// List of CommandGroup's
  std::vector<CommandGroup*> m_groups;

  /// Flag to know if method is to be called by a plugin and not by the subsystem
  bool m_isNonRootMethod;

public: // nested classes

  /// Class to allow STL algorithms to call member functions on this pointer
  /// This is similar to mem_fun_t from the STL standard
  template < typename Ret, typename TYPE >
    class root_mem_fun_t
    {
    public: // typedefs
      /// @c argument_type is the type of the argument
      typedef TYPE argument_type;
      /// @c result_type is the return type
      typedef Ret result_type;

    public: // functions
      /// constructor
      explicit root_mem_fun_t( Ret ( TYPE::*pf )() ) : m_f(pf) {}

      /// this is the operator we need to add in order for STL algorithms to function
      Ret operator()(TYPE* p) const
      {
        if(!p->isNonRootMethod()) return (p->*m_f)();
        return;
      }

    private: // data
      /// storage of the pointer to the member function to be called
      Ret (TYPE::*m_f)();

    }; // class root_mem_fun_t

}; // class Method

//////////////////////////////////////////////////////////////////////////////

/// Helper function to allow easy use of the root_mem_fun_t class
/// This is similar to mem_fun from the STL standard
template < typename Ret,  typename TYPE >
  inline typename Method::root_mem_fun_t<Ret,TYPE>
  root_mem_fun( Ret (TYPE::*f)() )
  {
    return typename Method::root_mem_fun_t<Ret,TYPE> (f);
  }

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Method_hh
