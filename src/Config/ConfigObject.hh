// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ConfigObject_hh
#define COOLFluiD_ConfigObject_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Common/NamedObject.hh"
#include "Common/FailedCastException.hh"

#include "Config/OptionList.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

/// Configurable object whose options can be set through a config file,
/// enviromental variables or by supplying command line arguments.
///
/// @author Tiago Quintino
class Config_API ConfigObject : public Common::NamedObject {

public: // functions

  /// Constructor
  ConfigObject(const std::string& name,
               const bool&   helpOnError = false,
               const std::string& usageSummary = std::string());

  /// Virtual destructor
  virtual  ~ConfigObject();

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure ( ConfigArgs& args );

  /// Configures the options for this object with a XML tree
  /// @param args is the XML tree to be parsed.
  virtual void configure_xml ( const std::string& xml_args );

  /// Configures the dynamic options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configureDynamic ( ConfigArgs& args );

  /// Configures the nested objects of this object with the arguments
  /// passed onto the object. This means that if used recursively all
  /// arguments are always passed and the configuration happens with
  /// both previuosly used and unused arguments.
  /// @param nestConfig is a pointer to the nested object to configure
  virtual void configureNested(ConfigObject* nestConfig, ConfigArgs& args );

  /// Configures the nested objects of this object with the arguments
  /// passed onto the object. This means that if used recursively all
  /// arguments are always passed and the configuration happens with
  /// both previuosly used and unused arguments.
  /// @param nestConfig is a pointer to the nested object to configure
  /// @param nestConfig is a reference to the nested object to configure
  virtual void configureNested(ConfigObject& nestConfig, ConfigArgs& args );

  /// Unconfigures the object
  /// Currently not in use.
  /// To be extended by derived classes.
  virtual void unconfigure();

  /// Check if the object is configured
  /// @return true if configure has been called
  virtual bool isConfigured() const { return m_config;}

  /// Provides a usage summary for all options of the object
  /// @return a string with the summary
  virtual std::string usageSummary() const;

  /// Returns the usage of this ConfigObject
  /// The format is that as for
  /// <pre>
  /// Usage: OBJECT NAME
  ///
  /// Options:
  /// name.opt     Option description (default = default)
  /// </pre>
  /// Note that there is one tab between the option and the description.
  /// @returns string with usage information
  std::string writeUsage() const;

  /// Prints processed and unprocessed arguments for debugging pruposes
  /// @returns string with debug information
  std::string debugInfo() const;

  /// Returns an XML with the tree of options and nested ConfigObject's
  /// @returns string with xml
  std::string getTreeXML() const;

  /// Returns whether there was an error during processing.
  /// @return true if there was an error
  virtual bool hadError() const;

  /// Sets the configuration file to be parsed for this object.
  /// The one added will follow any previously added configuration files.
  /// Note that the usage of "~" as the first character of the string will be converted
  /// to the value of the HOME environment variable.
  virtual void setConfigFile ( const boost::filesystem::path& file );

  /// Sets this object _nest member.
  /// @param nest the name of the parent object that nests this one
  void setNest(const std::string& nest);

  /// Gets this object _nest member.
  /// @return string with nest label
  std::string getNest() const;

  /// Gets this object full nested name.
  /// @return _nest + "." + _objName
  std::string getNestName() const;

  /// Sets the config options by calling the defineConfigOptions
  /// This will add nested names to the options as opposed to addOptionsTo
  /// @param prt should be passed with the this pointer to help identify callee the CLASS type
  template <typename CLASS>
  void addConfigOptionsTo(const CLASS* ptr)
  {
    CLASS::defineConfigOptions(m_options);
    registToConfigFacility();
  }

  /// Gets an option by its name
  /// @param name option name
  /// @return pointer to Option
  Option * getOption(const std::string& name) const;

  /// Gets an option by its name, but on its derived type
  /// @throw Common::FailedCastException if the cast to derived option fails
  /// @param name option name
  /// @return pointer to option derived type
  template < typename OPTIONTYPE >
  OPTIONTYPE * getOptionT(const std::string& name) const
  {
    Option * opt = getOption(name);
    OPTIONTYPE * dptr = dynamic_cast<OPTIONTYPE*>(opt);
    if (dptr == CFNULL)
    {
      std::string msg ("Option failed dynamic cast to ");
      msg += DEMANGLED_TYPEID(OPTIONTYPE);
      throw Common::FailedCastException (FromHere(),msg);
    }
    return dptr;
  }

  /// Associates the specified parameter to an option identified by its name.
  /// @param name option name
  /// @param value option value passed by pointer
  template < typename TYPE >
  void setParameter(const std::string& name, TYPE* const value)
  {
    Option * opt = getOption(name);
    // links the value to the option
    opt->linkToValue(value, DEMANGLED_TYPEID(TYPE));
    // sets the default value
    opt->setDefaultValue(value, DEMANGLED_TYPEID(TYPE));
  }

  /// Associates the specified parameter to an option identified by its name.
  /// Allows to set an option whose value type is different than the option type.
  /// @param name option name
  /// @param value option value passed by pointer
  /// @param default default value
  template < typename TYPE ,  typename DEFTYPE  >
  void setParameter(const std::string& name, TYPE* const value, DEFTYPE defvalue)
  {
    Option * opt = getOption(name);
    // links the value to the option
    opt->linkToValue(value, DEMANGLED_TYPEID(TYPE));
    // sets the default value
    opt->setDefaultValue(&defvalue, DEMANGLED_TYPEID(DEFTYPE));
  }

  /// Adds a nested ConfigObject
  void addNestedConfigObject(ConfigObject& obj);

  /// Dumps the Tree of options and Nested ConfigObject's to a string in XML format
  std::string getTreeXML();

  /// @returns the OptionList of this ConfigObject
  OptionList getOptionList() const { return m_options; }

  /// @returns the OptionList of this ConfigObject
  OptionList& getOptionList() { return m_options; }

  /// @returns the environmental variables
  std::vector<std::string> getEnvVars() const { return m_env_vars; }

  /// @returns the config file name
  boost::filesystem::path getConfigFile() const { return m_config_file; }

protected: // functions

  /// Processes the argument list from the configuration arguments
  /// Derived classes may override to change order in which options are parsed
  /// Default is:
  ///   1) config file via processConfigFile()
  ///   2) parent arguments via processParentArgs()
  ///   3) environmental variables via processEnvironmentVariables()
  /// Last configuration has preference
  virtual void processOptions ( ConfigArgs& args );

  /// Processes the argument list passed by client code
  /// Derived classes may change behavior, typically to skip altogether this configuration type
  virtual void processParentArgs ( ConfigArgs& args );

  /// Parses and processes the configuration file
  /// Derived classes may change behavior, typically to skip altogether this configuration type
  virtual void processConfigFile ( ConfigArgs& args );

  /// Processes the options defined by the environmental variables
  /// Derived classes may change behavior, typically to skip altogether this configuration type
  virtual void processEnvironmentVariables ( ConfigArgs& args );

  /// Reads the environment variable, splits it by whitespace,
  /// and processes it as a command line.
  virtual std::pair<bool,std::string> readEnvironmentVariable(const std::string& envVar);

  /// Adds an environment variable to be parsed.
  /// Will follow any previously added environment variables
  virtual void addEnvironmentVariable(const std::string& envVar);

  /// Regists this object to the ConfigFacility
  void registToConfigFacility();

  /// Removes this object registration from the ConfigFacility
  void unregistFromConfigFacility();

private: // member data

  /// the list of options of this object
  OptionList m_options;

  /// the list of nested config objects
  std::vector< ConfigObject* > m_nested_objs;

  /// Setup flag to indicate that the object is configurated.
  bool m_config;

  /// Whether to show help when processing fails.
  bool m_helpOnError;

  /// Whether there was an error when processing the arguments.
  bool m_hadError;

  /// Whether the object has been registered to the ConfigFacility
  bool m_registered;

  /// Configuration file
  boost::filesystem::path m_config_file;

  /// The environment variables.
  std::vector<std::string> m_env_vars;

  /// The nested names of this object.
  std::string m_nest;

  /// Summary of usage.
  std::string m_usageSummary;

  /// Pointer to the parent ConfigObject
  ConfigObject * parentConfig;

 private: // functions

  /// Sets the ConfigObject given by parameter as the parent of this ConfigObject
  /// object. If this object already has a parent, nothing is done.
  /// \param parentConfig Pointer to the parent ConfigObject
  void setParentConfig(ConfigObject * parentConfig);

  /// Unregisters (erases from the vector) the given child ConfigObject object.
  /// If the child does not belong to this object, nothing is done.
  /// \param config The child to unregister.
  void unregisterFromParent(ConfigObject * config);

public: // member data

  /// Separator for nested configs
  static const std::string NEST_SEPARATOR;

  /// Separator for nested configs but for the enviromental variables
  /// because there are some characters that are not allowed, e.g. the dot.
  static const std::string NEST_SEPARATOR_ENVVAR;

}; // class ConfigObject

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ConfigObject_hh
