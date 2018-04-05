// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MethodData_hh
#define COOLFluiD_Framework_MethodData_hh

#include "Common/SafePtr.hh"
#include "Common/OwnedObject.hh"
#include "Common/SetupObject.hh"
#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"
#include "Framework/NamespaceMember.hh"
#include "Framework/CollaboratorAccess.hh"
#include "Framework/MethodStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Method;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the Data aggregator of a method
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API MethodData : public Common::NonCopyable<MethodData>,
                   public Common::OwnedObject,
                   public Common::SetupObject,
                   public Config::ConfigObject,
                   public Framework::NamespaceMember,
                   public Framework::CollaboratorAccess {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures this MethodData
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Returns the list of strategies of tyhis Method
  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategies() const 
  { return m_strategies; }
  
  /// Template function for creating and configuring MethodStrategy's with reduced
  /// checking of template parameters because BASECOMMAND can be more specific than a MethodStrategy.
  /// @param stg  the pointer which should be holding the newly created and configured strategy
  /// @param type the string with the name of the type of strategy to create
  /// @param name the name to give to the newly created strategy
  /// @param data the MethodData that the strategy will share with other commands and strategies of the same method
  template < typename BASESTRATEGY, typename DATA>
  void configureStrategy( Config::ConfigArgs& args,
                          Common::SelfRegistPtr<BASESTRATEGY>& stg,
                          const std::string& type,
                          const std::string& name,
                          const Common::SharedPtr<DATA>& data )
  {
    typedef typename BASESTRATEGY::PROVIDER ProviderType;
    CFLogDebugMed( "MethodData::configureStrategy() => Type: " << type << " Name: " << name << "\n");
    Common::SafePtr<ProviderType> prov;
    try
    {
      prov = FACTORY_T_GET_PROVIDER(this->getFactoryRegistry(), BASESTRATEGY, type);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      throw;
    }
    cf_assert(prov.isNotNull());
    stg = prov->create(name,data);
    cf_assert(stg.isNotNull());
    stg->setFactoryRegistry(this->getFactoryRegistry());
    m_strategies.push_back(stg.getPtr());
    configureNested ( stg.getPtr(), args );
  }

  /// Template function for creating and configuring NumericalStrategy's
  /// No MethodData is required
  /// @param stg  the pointer which should be holding the newly created and configured strategy
  /// @param type the string with the name of the type of strategy to create
  /// @param name the name to give to the newly created strategy
  /// @param data the MethodData that the strategy will share with other commands and strategies of the same method
  template < typename BASESTRATEGY >
  void configureStrategy ( Config::ConfigArgs& args,
                           const std::string& type,
                           const std::string& name,
                           Common::SelfRegistPtr<BASESTRATEGY>& stg )
  {
    typedef typename BASESTRATEGY::PROVIDER ProviderType;
    CFLogDebugMed( "MethodData::configureStrategy() => Type: " << type << " Name: " << name << "\n");
    Common::SafePtr<ProviderType> prov;
    try
    {
      prov = FACTORY_T_GET_PROVIDER(this->getFactoryRegistry(), BASESTRATEGY, type);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      throw;
    }
    cf_assert(prov.isNotNull());
    stg = prov->create(name);
    cf_assert(stg.isNotNull());
    stg->setFactoryRegistry(this->getFactoryRegistry());
    m_strategies.push_back(stg.getPtr());
    configureNested ( stg.getPtr(), args );
  }

  /// Template function for creating and configuring MethodStrategy's with reduced
  /// checking of template parameters because BASECOMMAND can be more specific than a MethodStrategy.
  /// @param stg  the pointer which should be holding the newly created and configured strategy
  /// @param type the string with the name of the type of strategy to create
  /// @param name the name to give to the newly created strategy
  /// @param data the MethodData that the strategy will share with other commands and strategies of the same method
  template < typename BASESTRATEGY, typename DATA >
  void configureStrategy( Config::ConfigArgs& args,
                          const std::string& type,
                          const std::string& name,
                          Common::SelfRegistPtr<BASESTRATEGY>& stg,
                          const Common::SharedPtr<DATA>& data )
  {
    typedef typename BASESTRATEGY::PROVIDER ProviderType;
    CFLogDebugMed( "MethodData::configureStrategy() => Type: " << type << " Name: " << name << "\n");
    Common::SafePtr<ProviderType> prov;
    try
    {
      prov = FACTORY_T_GET_PROVIDER(this->getFactoryRegistry(), BASESTRATEGY, type);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      throw;
    }
    cf_assert(prov.isNotNull());
    stg = prov->create(name);
    cf_assert(stg.isNotNull());
    stg->setFactoryRegistry(this->getFactoryRegistry());
    stg->setMethodData(data);
    m_strategies.push_back(stg.getPtr());
    configureNested ( stg.getPtr(), args );
  }

  /// Template function for creating and configuring MethodStrategy's with reduced
  /// checking of template parameters because BASECOMMAND can be more specific than a MethodStrategy.
  /// @param stg  the pointer which should be holding the newly created and configured strategy
  /// @param type the string with the name of the type of strategy to create
  /// @param name the name to give to the newly created strategy
  /// @param data the MethodData that the strategy will share with other commands and strategies of the same method
  template < typename BASESTRATEGY, typename DATA, typename ARG1 >
  void configureStrategy( Config::ConfigArgs& args,
                          const std::string& type,
                          const std::string& name,
                          Common::SelfRegistPtr<BASESTRATEGY>& stg,
                          const Common::SharedPtr<DATA>& data,
                          const ARG1& arg1 )
  {
    typedef typename BASESTRATEGY::PROVIDER ProviderType;
    CFLogDebugMed( "MethodData::configureStrategy() => Type: " << type << " Name: " << name << "\n");
    Common::SafePtr<ProviderType> prov;
    try
    {
      prov = FACTORY_T_GET_PROVIDER(this->getFactoryRegistry(), BASESTRATEGY, type);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      throw;
    }
    cf_assert(prov.isNotNull());
    stg = prov->create(name,arg1);
    cf_assert(stg.isNotNull());
    stg->setFactoryRegistry(this->getFactoryRegistry());
    stg->setMethodData(data);
    m_strategies.push_back(stg.getPtr());
    configureNested ( stg.getPtr(), args );
  }

protected: // functions
  
  /// Get the provider corresponding to the given name
  template <typename OBJ>
  Common::SafePtr<typename OBJ::PROVIDER> getProvider(std::string& provName)
  {
    Common::SafePtr<typename OBJ::PROVIDER> prov = CFNULL;
    try {
      prov = FACTORY_T_GET_PROVIDER(this->getFactoryRegistry(), OBJ, provName);
      CFLog(VERBOSE, this->getName() << "::configure() => Provider " << provName << " exists\n");
    }
    catch (Common::NoSuchValueException& e) {
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing Null for " << OBJ::getClassName() << " instead ..." << "\n");
      prov = FACTORY_T_GET_PROVIDER(getFactoryRegistry(), OBJ, "Null");
    }
    return prov;
  }    
  
  /// Default constructor without arguments.
  MethodData(Common::SafePtr<Method> owner);

  /// Destructor.
  virtual ~MethodData();

protected: // data

  /// storage of the strategies in this method
  /// @todo transfer ownership of the Strategies here
  std::vector <Common::SafePtr<NumericalStrategy> > m_strategies;

}; // end of class MethodData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MethodData_hh
