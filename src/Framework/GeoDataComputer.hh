// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GeoDataComputer_hh
#define COOLFluiD_Framework_GeoDataComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/Node.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for (re-)computing algorithm-dependent
/// geometric data (normals, volumes, etc.), assuming that all data arrays have
/// been already resized correctly
/// @author Andrea Lani
template < typename METHODDATA >
class GeoDataComputer : public Framework::MethodStrategy<METHODDATA> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider<METHODDATA,GeoDataComputer<METHODDATA> > PROVIDER;

  /// Constructor
  GeoDataComputer(const std::string& name);

  /// Default destructor
  virtual ~GeoDataComputer();

  /// Set private data that will be used during the computation
  virtual void setup() 
  {
    Framework::MethodStrategy<METHODDATA>::setup();
  }
  
  /// Compute the geometric data and fill pre-allocated data arrays 
  virtual void compute() = 0;
   
  /// Unsetup private data
  virtual void unsetup()
  {
    Framework::MethodStrategy<METHODDATA>::unsetup();
  }
  
  /// Compute the nodes at intermediate time for the cell centers
  virtual void modifyOffMeshNodes() = 0;
  
  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }
  
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<BaseDataSocketSink> > result;
    result.push_back(&socket_states);
    result.push_back(&socket_nodes);
    return result;
  }
  
  /// Gets the Class name
  static std::string getClassName() {return "GeoDataComputer";}

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

protected:

  /// Private Copy Constructor
  GeoDataComputer(const GeoDataComputer& o);

  /// Private Assignement operator
  const GeoDataComputer& operator=(const GeoDataComputer& o);
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
}; // end of class GeoDataComputer
      
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "GeoDataComputer.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GeoDataComputer_hh
