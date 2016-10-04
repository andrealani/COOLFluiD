// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataHandleOutput_hh
#define COOLFluiD_Framework_DataHandleOutput_hh

//////////////////////////////////////////////////////////////////////////////

#include <ostream>

#include "Common/Trio.hh"
#include "Common/NonCopyable.hh"
#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/NumericalStrategy.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class outputs data from a datahandle of reals
/// @todo could be made a template and provide type independent output
/// @author Tiago Quintino
class Framework_API DataHandleOutput :
  public NumericalStrategy,
  public Framework::NamespaceMember {

public: // typedefs

  typedef Environment::ConcreteProvider<DataHandleOutput,1> PROVIDER;
  typedef const std::string& ARG1;

  typedef Common::Trio < CFuint, CFuint, DataHandle<CFreal> > DataHandleInfo;

public: // functions

  /// Constructor
  DataHandleOutput(const std::string& name);

  /// Virtual destructor
  virtual ~DataHandleOutput();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup the object
  virtual void setup();

  /// Unsetup this object
  virtual void unsetup();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "DataHandleOutput"; }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
    
  /// Outputs the names of the variables
  virtual void printVarNames (std::ostream& out) const ;

  /// prints the data which is assumed to be relative to id of a state
  /// @pre must call first getDataHandles
  virtual void printStateData (std::ostream& out, CFuint state_id) const ;
  
  /// fill a given array with the data handle data corresponding to a given state ID
  /// @post uses and updates the given counter for indexing the output array (out)
  virtual void fillStateData(CFreal* out, CFuint state_id, CFuint& counter, int iVar = -1) const;
  
  /// prints the data which is assumed to be relative to id of a cell
  /// @pre must call first getDataHandles
  virtual void printCCData (std::ostream& out, CFuint cell_id) const ;
 
  /// fill a given array with the data handle data corresponding to a given state ID
  /// @post uses and updates the given counter for indexing the output array (out)
  virtual void fillStateDataCC(CFreal* out, CFuint cell_id, CFuint& counter, int iVar = -1) const;
  
  /// gets the raw data corresponding to the variable index of state based variables
  /// @param var_id should be less than m_varnames.size()
  DataHandleInfo getCCData(CFuint var_id) const;

  /// gets the raw data for the variable index
  /// @param var_id should be less than m_varnames.size()
  DataHandleInfo getStateData(CFuint var_id) const;

  /// prints the data which is assumed to be relative to id of a state
  virtual void getDataHandles ();

  /// gets the names for each variable
  virtual std::vector<std::string> getVarNames () const;

  /// gets the names for each variable for the cell centered datahandles
  virtual std::vector<std::string> getCCVarNames () const;

  /// gets the names of the Trs which hold the elements for CC datahandles
  virtual std::vector<std::string> getCCVarTrs () const;

  /// Returns the DataSocket's that this strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// helper function to setup the variables based on the states
  virtual void setupStateVars();

  /// helper function to setup the variables based on the cell center
  virtual void setupCCVars();

private: // data

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// the dynamic sockets to access the data to output
  Framework::DynamicDataSocketSet<> m_sockets;

  /// socket names ofr state variables
  std::vector<std::string> m_socket_names;

  /// socket names for cell-centered values
  std::vector<std::string> m_ccsocket_names;

  /// nb of state variables per socket
  std::vector<CFuint> m_nbvars;

  /// nb variables per socket for cell-centered values
  std::vector<CFuint> m_ccnbvars;

  /// name of each state variable
  std::vector<std::string> m_varnames;

  /// name of each cell-centered variable
  std::vector<std::string> m_ccvarnames;

  /// name of the TRS from where to output the cells
  /// data will be zero for the other TRS's
  /// @pre the TRS must be tagged "inner" and "cells"
  /// @post sized with number of cc variables
  std::vector<std::string> m_ccvartrs;

  /// name of the TRS from where to output the cells
  /// data will be zero for the other TRS's
  /// @pre the TRS must be tagged "inner" and "cells"
  /// @post sized with number of cc sockets
  std::vector<std::string> m_ccsockettrs;

  /// cache of datahandles, but gets refreshed every iteration
  std::vector< DataHandle<CFreal> > m_datahandles;

  /// cache of datahandles, but gets refreshed every iteration
  std::vector< DataHandle<CFreal> > m_ccdatahandles;

  /// is the cache of datahandles ready ?
  bool m_cached_datahandles;

  /// block size for each cell-centered variable
  /// this is usefull to circunvent the existence of subcells
  /// in that case set this value to the maximum number of subcells
  std::vector< CFuint > m_ccblocksize;
  
  /// flag telling if the variables are nodal
  bool m_isNodal;
  
}; // end of class DataHandleOutput

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(DataHandleOutput) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataHandleOutput_hh
