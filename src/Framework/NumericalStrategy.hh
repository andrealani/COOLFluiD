// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NumericalStrategy_hh
#define COOLFluiD_Framework_NumericalStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/SafePtr.hh"
#include "Common/OwnedObject.hh"
#include "Common/SetupObject.hh"
#include "Config/ConfigObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class BaseDataSocketSource;
    class BaseDataSocketSink;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalStrategy.
/// It is the abstract class definning the interface to the derived NumericalStrategys
/// NumericalStrategy allow a given Method to be costumized through the setup
/// of the commands without subclassing. A given Method might use a several ways to compute
/// something. Instead of multiclassing, the user chooses a suitable Strategy and configures
/// the method with it.
/// @author Tiago Quintino
class Framework_API NumericalStrategy : public Common::OwnedObject,
                          public Common::SetupObject,
                          public Config::ConfigObject,
                          public Common::NonCopyable<NumericalStrategy>     {
public:

  typedef Environment::ConcreteProvider<NumericalStrategy> PROVIDER;

  /// Default constructor
  explicit NumericalStrategy(const std::string& name);

  /// Default destructor pure virtual
  virtual ~NumericalStrategy() = 0;

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unset up private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup();

  /// Configure
  virtual void configure ( Config::ConfigArgs& args );

  /// Configure the nested sockets in this NumericalStrategy
  void configureNestedSockets ( Config::ConfigArgs& args );

  /// Allocates all sockets in this NumericalStrategy
  void allocateStrategySockets();

  /// Deallocates all sockets in this NumericalStrategy
  void deallocateStrategySockets();

  /// Returns the DataSocket's that this strategy provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<BaseDataSocketSink> > needsSockets();

  /// Gets the Class name
  static std::string getClassName() { return "NumericalStrategy"; }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() = 0;
  
}; // class NumericalStrategy

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(NumericalStrategy)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NumericalStrategy_hh
