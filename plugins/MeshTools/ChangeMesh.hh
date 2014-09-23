// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ChangeMesh_hh
#define COOLFluiD_MeshTools_ChangeMesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "MeshTools/QualityCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////
		
struct MeshChangeSet 
{
  
  // Data defining the mesh change to pass to algorithm
  
  /// Indexes of nodes to remove (ordered)
	std::vector < CFuint > remove_nodes; 
  /// Indexes of states to remove (ordered)
	std::vector < CFuint > remove_states;
	
	
  /// global ID of state to be added (ordered)
	std::vector < CFuint > add_states_gid;
  /// Coordinates of states to be added
	std::vector < RealVector > add_states_coords;
  /// Values of states to be added
	std::vector < RealVector > add_states_values;
  /// if the state to be added is an updatable state
	std::vector < bool > add_states_updatable;
  
  
  /// global ID of nodes to be added (ordered)
	std::vector < CFuint > add_nodes_gid;
  /// Coordinates of nodes to be added
	std::vector < RealVector > add_nodes_coords;
  /// if the node to be added is an updatable state
	std::vector < bool > add_nodes_updatable;
  
  
  // Data return from mesh change algorithm

  /// returns the idx of added nodes here (ordered)
  /// will be resized by algorithm
  std::vector < CFuint > added_nodes_local_idx;
  /// returns the idx of added states here (ordered)
  /// will be resized by algorithm
  std::vector < CFuint > added_states_local_idx;
	
};

//////////////////////////////////////////////////////////////////////////////

/// This class computes changes the mesh and afterwards 
/// raises an event to all methods resetup
/// @author Tiago Quintino
class ChangeMesh : public Framework::DataProcessingCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  ChangeMesh(const std::string& name);

  /// Virtual destructor
  virtual ~ChangeMesh();

  /// Setup private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unsetup private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void unsetup();

  /// Executes this command
  virtual void execute();

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  void add_remove_nodes ( MeshChangeSet& chset );

  void add_remove_states ( MeshChangeSet& chset );

  void add_remove_geoents ( MeshChangeSet& chset );

private: //data

  /// storage of the States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of the States
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

}; // end of class ChangeMesh

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ChangeMesh_hh
