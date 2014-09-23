// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_OptionMethodStrategy_hh
#define COOLFluiD_Framework_OptionMethodStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/OptionList.hh"
#include "Framework/BaseMethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

///  Handles configration of method strategy
///  @author Tiago Quintino
template < typename MSTRAT, typename BASETYPE >
class OptionMethodStrategy : public Config::Option {

public: // functions

  /// Constructor
  OptionMethodStrategy(const std::string& name, const std::string& description) :
    Config::Option(name, description), m_default_str(),  m_config_str(), m_ptr(CFNULL)  {}

  /// Copy constructor
  OptionMethodStrategy(const OptionMethodStrategy<MSTRAT,BASETYPE>& rhs) : Config::Option(rhs)
  {
    operator=(rhs);
  }

  /// Assignement operator
  OptionMethodStrategy<MSTRAT,BASETYPE>& operator=(const OptionMethodStrategy<MSTRAT,BASETYPE>& rhs)
  {
    Option::operator=(rhs);
    m_default_str = rhs.m_default_str;
    m_config_str  = rhs.m_config_str;
    m_ptr         = rhs.m_ptr;
    return *this;
  }

  /// Virtual destructor.
  virtual ~OptionMethodStrategy() {}

  void putMethodData(MSTRAT* data) { m_mdata = data; }

  void linkToValue(void* value, const std::string& typestr)
  {
    cf_assert( value != CFNULL );

    std::string thistype = DEMANGLED_TYPEID(Common::SelfRegistPtr<BASETYPE>);
    if( thistype != typestr )
      throw Config::BadMatchException
      (FromHere(), "A type [" + typestr + "] was passed to option [" +
        getLabel()  + "] which requires a value of type [" + thistype + "]" );

    m_ptr = static_cast<Common::SelfRegistPtr<BASETYPE>*>(value);
  }

  void setDefaultValue(void* defvalue, const std::string& typestr)
  {
    cf_assert( defvalue != CFNULL );

    std::string thistype = DEMANGLED_TYPEID(std::string);
    if( thistype != typestr )
      throw Config::BadMatchException
      (FromHere(), "A type [" + typestr + "] was passed as default to option [" +
        getLabel()  + "] which requires a value of type [" + thistype + "]" );

    std::string * ptr_value  = static_cast<std::string*>(defvalue);
    m_default_str = *ptr_value;
  }

  /// Sets this option value from a std::string
  /// @pre label has already been checked to ensure that it matches the value's tag
  /// @param value string with value to be converted to actual type
  virtual void setValue(const std::string& value)
  {
    CFLogDebugMed("Setting method strategy value [" << value << "]\n");

    m_config_str = value;
    this->validate();
    setLinkedOptions(value);
  }

  /// Returns the default value, as a std::string
  virtual std::string getDefaultValueAsString() const
  {
    return m_default_str;
  }

  /// Returns the default value, as a std::string
  virtual std::string getValueAsString() const
  {
    return m_config_str;
  }

  /// Returns the low-level details about this object
  /// @return string with details
  virtual std::string debugInfo() const
  {
    return "{ " + getLabel() + ", " + getHierarchyOwner() +  ", " + description() + m_config_str + " }";
  }

  /// Returns options in XML format
  /// @return string with description
  virtual std::string getXML() const
  {
      std::ostringstream os;
      os << "<" << getLabel() << "\t";
      os << "tree=\"" << "object" << "\" ";
      os << "type=\"" << DEMANGLED_TYPEID(MSTRAT) << "\" ";
      os << "abstype=\"" << DEMANGLED_TYPEID(BASETYPE) << "\" ";
      os << "mode=\"" << ( isBasic() ? "basic" : "advanced" ) << "\" ";
      os << "dynamic=\"" << ( isDynamic() ? "dynamic" : "static" ) << "\" ";
      os << "description=\"" << description() << "\" ";
      os << ">\n";
      os << getValueAsString() << "\n";
      os << "</" << getLabel() << ">\n";
      return os.str();
  }

  /// Sets values from command line parameters
  virtual CFuint setCommandLineValues(const std::string& arg,
                                      CFuint position,
                                      int argc,
                                      char** argv)
  {
      int nArgs = 1; ++position;
      if (position < CFuint(argc))
      {   // there needs to be a following argument, which is the value
          std::string argstr(argv[position]);
          checkArgument(argstr);
          setValue(argstr);
          position++;         // and now consume the value
          nArgs++;
      }
      else { throw Config::BadMatchException (FromHere(),"argument expected for " + getLabel()); }
      return nArgs;
  }

  /// Creates a new Option* with the same values as this
  virtual Option* clone() const { return new OptionMethodStrategy<MSTRAT,BASETYPE>(*this); }

  /// Called after all options have been successfully processed.
  virtual void onComplete(const Config::OptionList& opts)
  {
    Option::onComplete(opts);

    std::string provname = m_config_str.empty() ? m_default_str : m_config_str ;

    cf_assert_desc("MethodData must put itself into the option", m_mdata != CFNULL);

    CFLogDebugMed("Configuring method strategy ["
                  << BASETYPE::getClassName()
                  << "] with  value ["
                  << provname << "]\n");

    Common::SafePtr<BaseMethodStrategyProvider<MSTRAT,BASETYPE> > prov =
      Environment::Factory<BASETYPE>::getInstance().getProvider(provname);
    cf_assert(prov.isNotNull());

    (*m_ptr) = prov->create(provname, Common::SharedPtr<MSTRAT>(m_mdata));
    cf_assert(m_ptr->isNotNull());
  }

  /// Called when the value is set interactivelly
  /// @return true if successful
  virtual void onUpdate()
  {
    Option::onUpdate();

    std::string old_prov_name = m_ptr->getProviderBase()->getProviderName();

    std::string provname = m_config_str.empty() ? m_default_str : m_config_str ;

    cf_assert_desc("MethodData must put itself into the option",m_mdata != CFNULL);
    cf_assert_desc("Provider name must be different when calling onUpdate", provname != old_prov_name);

    CFLogWarn("Updating the method strategy [" << BASETYPE::getClassName()
              << "] from  value [" << old_prov_name << "]"
              << "  to value [" << provname << "]\n");

    Common::SafePtr<BaseMethodStrategyProvider<MSTRAT,BASETYPE> > prov =
      Environment::Factory<BASETYPE>::getInstance().getProvider(provname);
    cf_assert(prov.isNotNull());

    (*m_ptr) = prov->create(provname, Common::SharedPtr<MSTRAT>(m_mdata));
    cf_assert(m_ptr->isNotNull());
  }

  /// Check if option has been set and configured
  bool isConfigured () const { return m_ptr->isNotNull(); }

private: // data

  /// The default value as string
  std::string m_default_str;

  /// The string for configuration
  std::string m_config_str;

  /// The pointer to be configured
  Common::SelfRegistPtr<BASETYPE> * m_ptr;

  /// the pointer to the method data
  MSTRAT * m_mdata;

}; // class OptionMethodStrategy

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_OptionMethodStrategy_hh
