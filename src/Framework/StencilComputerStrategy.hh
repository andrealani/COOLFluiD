// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StencilComputerStrategy_hh
#define COOLFluiD_Framework_StencilComputerStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodStrategy.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class provides an abstract interface for computing
   * stencils
   * 
   * @author Willem Deconinck
   */
template <typename METHODDATA>
class StencilComputerStrategy : public Framework::MethodStrategy<METHODDATA> {

public:

  typedef Framework::BaseMethodStrategyProvider
  <METHODDATA,StencilComputerStrategy<METHODDATA> > PROVIDER;

  /**
   * Defines the config options of this class
   * @param   options   config options of this class
   */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  StencilComputerStrategy(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~StencilComputerStrategy();

  /** 
   * Compute the stencil
   */
  virtual void compute () =0;
  virtual void computeInspected () =0;
  virtual void computeStencil (const CFuint& iStencil) =0;
  
  /**
   * Output Stencil
   */
   virtual void outputStencil () =0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "StencilComputerStrategy";
  }

  /** 
   * Set private data that will be used during the computation
   */
  virtual void setup();

  /** 
   * Configure the object
   * @param    args   configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args );

  /** 
   * Returns the DataSockets that this strategy needs 
   * as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
    needsSockets();

  /**
   * Returns the DataSockets that this strategy provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
    providesSockets();

protected: //data

  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

}; // end of class ComputeStencil

//////////////////////////////////////////////////////////////////////////////

template <typename METHODDATA>
void StencilComputerStrategy<METHODDATA>::defineConfigOptions(Config::OptionList& options)
{
  // Configuration option definitions here	
}

//////////////////////////////////////////////////////////////////////////////

template <typename METHODDATA>
StencilComputerStrategy<METHODDATA>::StencilComputerStrategy(const std::string& name) :
  Framework::MethodStrategy<METHODDATA>(name),
    socket_states("states"),
    socket_nodes("nodes")
{
  this->addConfigOptionsTo(this);
  // Setting default configurations here.
}

//////////////////////////////////////////////////////////////////////////////

template <typename METHODDATA>
StencilComputerStrategy<METHODDATA>::~StencilComputerStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

template <typename METHODDATA>
void StencilComputerStrategy<METHODDATA>::setup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

template <typename METHODDATA>
void StencilComputerStrategy<METHODDATA>::configure ( Config::ConfigArgs& args )
{
  Framework::MethodStrategy<METHODDATA>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

template <typename METHODDATA>
std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > StencilComputerStrategy<METHODDATA>::needsSockets()
{
  // create socket sink
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

template <typename METHODDATA>
std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > StencilComputerStrategy<METHODDATA>::providesSockets()
{
  // create socket source
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_StencilComputerStrategy_hh
