// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_EquationFilter_hh
#define COOLFluiD_Framework_EquationFilter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for an equation filter
/// @author Andrea Lani
template < typename METHODDATA >
class EquationFilter : public Framework::MethodStrategy<METHODDATA> {
public:
  
  typedef Framework::BaseMethodStrategyProvider
  <METHODDATA,EquationFilter<METHODDATA> > PROVIDER;
  
  /// Constructor
  EquationFilter(const std::string& name) :
    MethodStrategy<METHODDATA>(name),
    socket_states("states"),
    socket_nodes("nodes"),
    socket_activeDiffusion("activeDiffusion")
  {
    this->addConfigOptionsTo(this);
    
    m_startIter = 0.;
    this->setParameter("StartIter",&m_startIter);
  }
  
  /// Default destructor
  virtual ~EquationFilter()
  {
  }
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options)
  {
    options.addConfigOption< CFuint, Config::DynamicOption<> >("StartIter","Iteration at which filtering must start.");
  }
  
  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }
  
  /// Set up private data to prepare the simulation
  virtual void setup() 
  {
    Framework::MethodStrategy<METHODDATA>::setup();
  }
  
  /// Unsetup private data to prepare the simulation
  virtual void unsetup() 
  {
    Framework::MethodStrategy<METHODDATA>::unsetup();
  }
  
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    using namespace std;
    
    std::vector<Common::SafePtr<BaseDataSocketSink> > result;
    result.push_back(&socket_states);
    result.push_back(&socket_nodes);
    result.push_back(&socket_activeDiffusion);
    return result;
  }
  
  /// Reset all data before starting a new iteration
  virtual void reset() = 0;
  
  /// Filter the equation on the current geometric entity
  virtual bool filterOnGeo(GeometricEntity *const geo) = 0;
  
  /// Gets the Class name
  static std::string getClassName() {  return "EquationFilter"; }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
private: // methods

  /// Private Copy Constructor
  EquationFilter(const EquationFilter& f);

  /// Private Assignement operator
  const EquationFilter& operator=(const EquationFilter& f);
  
protected:
  
  /// storage of the states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// storage of the nodes
  Framework::DataSocketSink<Framework::Node*,  Framework::GLOBAL> socket_nodes;
  
  /// flags telling if diffusion terms have to be computed for a given cell
  Framework::DataSocketSink<CFreal> socket_activeDiffusion;
  
  /// Iteration at which filtering must start
  CFuint m_startIter;
 
}; // end of class EquationFilter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_EquationFilter_hh
