// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/EventHandler.hh"
#include "Environment/DirPaths.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "MeshTools/MeshTools.hh"
#include "MeshTools/ChangeMesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ChangeMesh, DataProcessingData, MeshToolsModule> ChangeMeshProvider("ChangeMesh");

//////////////////////////////////////////////////////////////////////////////

void ChangeMesh::defineConfigOptions(Config::OptionList& options)
{
//   options.addConfigOption< std::string >("QualityType","Type of quality to use for the measure");
}

//////////////////////////////////////////////////////////////////////////////

ChangeMesh::ChangeMesh(const std::string& name) : 
  DataProcessingCom(name),
  socket_states("states"),
  socket_nodes("nodes")
{
  CFAUTOTRACE;

	addConfigOptionsTo(this);

//  m_qualityType = "Concrete";
//	setParameter("QualityType",&_qualityType);
  
  CFLog ( NOTICE, "++++++ ChangeMesh::ChangeMesh() [" << getName() << "]\n" );
}

//////////////////////////////////////////////////////////////////////////////

ChangeMesh::~ChangeMesh()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ChangeMesh::setup()
{
  CFAUTOTRACE;
  
  DataProcessingCom::setup();
  
  CFLog ( NOTICE, "++++++ ChangeMesh::setup() [" << getName() << "]\n" );
}

//////////////////////////////////////////////////////////////////////////////

void ChangeMesh::unsetup()
{
  CFAUTOTRACE;
  
  CFLog ( NOTICE, "++++++ ChangeMesh::unsetup() [" << getName() << "]\n" );
  
  DataProcessingCom::unsetup();
}
    
//////////////////////////////////////////////////////////////////////////////

void ChangeMesh::add_remove_nodes( MeshChangeSet& chset )
{
  CFAUTOTRACE;
  
  DataHandle < Node*, GLOBAL >  nodes = socket_nodes.getDataHandle();
  
  
  // some initial sanity checks
  
  const CFuint nb_nodes_to_remove = chset.remove_nodes.size();
  cf_assert ( chset.add_nodes_coords.size()    == chset.add_nodes_gid.size() );
  cf_assert ( chset.add_nodes_updatable.size() == chset.add_nodes_gid.size() );
  
  
  /// @TODO note that we are not deleting the nodes that are removed
  ///       we reuse their memory instead of adding a new one, but those 
  ///       that are left are not deleted
  
  const CFuint nb_nodes_to_increase = 
  chset.add_nodes_gid.size() > chset.remove_nodes.size() ? 
  chset.add_nodes_gid.size() - chset.remove_nodes.size() : 0 ;
  
  // remember the old number of nodes
  const CFuint old_nodes_size = nodes.size();
  CFLog ( NOTICE , "OLD Nb Nodes [" << old_nodes_size << "]\n" );
  
    
  
  // resize the container of nodes
  nodes.resize( old_nodes_size + nb_nodes_to_increase );
  const CFuint new_nodes_size = nodes.size();
  
  
  // correctly size the indexes of the nodes to add
  chset.added_nodes_local_idx.resize(chset.add_nodes_gid.size());
  
    
  
  // loop over the old nodes to remove nodes
  // and if needed add new ones in their place reusing their local ID
  
  CFuint last_removed_ndix = 0; // will vary from 0 to nb nodes to remove
  CFuint last_added_ndix = 0;   // will vary from 0 to nb nodes to add
  
  for ( CFuint i = 0;  i < old_nodes_size; ++i )
  {
      
    
    CFLog ( NOTICE, "NODE [" << i << "]\n");
    
    bool node_removed = false;
    
      
    
    // verify that we are not passing ahead of node indexes to remove
    // we assume that remove_nodes is an ordered vector
    if ( nb_nodes_to_remove != 0 )
    {    
      cf_assert ( i <= chset.remove_nodes[last_removed_ndix] );
    }
    
    
    // is this node to be removed?
    if ( last_removed_ndix < nb_nodes_to_remove && chset.remove_nodes[last_removed_ndix] == i )
    {
      
  
      
      // node memory is not removed
      node_removed = true;
      
      CFLog ( NOTICE, "Remove node with local index [" << i << "]\n");
      
      // incremet index of nodes to remove
      ++last_removed_ndix; 
    }
    
      
    
    // reuse the removed index if there are nodes to add
    if ( node_removed && last_added_ndix < chset.add_nodes_gid.size() ) 
    {
  
      
      const CFuint node_global_id = chset.add_nodes_gid[last_added_ndix];
      const CFuint node_local_id  = i;
      
      chset.added_nodes_local_idx[last_added_ndix] = node_local_id;
      
      
      CFLog ( NOTICE, "Replaced with node with global index [" << node_global_id << "]\n");
      
      const bool isUpdatable = chset.add_nodes_updatable[last_added_ndix];
      
      nodes[i]->setParUpdatable( isUpdatable );
      nodes[i]->setGlobalID( node_global_id );
      nodes[i]->setLocalID ( node_global_id );
      
      // put the values in the added node
      // verify the values are correctly sized
      RealVector& new_node_coords = chset.add_nodes_coords[last_added_ndix];
      cf_assert ( nodes[i]->size() == new_node_coords.size() );
      
      *nodes[i] = new_node_coords;
      
      // increment index of added nodes
      ++last_added_ndix;     
    }
    
      
    
  }
  
    
  
  // check all nodes where removed
  cf_assert ( last_removed_ndix == chset.remove_nodes.size() );
  
  // check all nodes where added
  // if not we need to add more points
  if ( last_added_ndix < chset.add_nodes_gid.size() )
  {
    
    
    
    // start from the last node index
    CFuint new_node_idx = old_nodes_size;
    
    // loop over the remaining nodes to add
    // and grow the nodes memory
    for ( ; last_added_ndix < chset.add_nodes_gid.size(); ++last_added_ndix )
    {
      // global id of node to add
      const CFuint node_global_id = chset.add_nodes_gid[last_added_ndix];
      
      CFLog ( NOTICE, "NEW NODE [" << last_added_ndix << "] with GLOBAL ID [" << node_global_id << "]\n");
      

      
      // get the node_local_id from the global ID
      const CFuint node_local_id  = nodes.addLocalPoint ( node_global_id );
      

      
      // check that the returned local ID matches 
      // the new node idx ccording to our accounting
      cf_assert ( new_node_idx == node_local_id );
      
      // we already resize the nodes to the new size so the index 
      // should be lower that the size
      cf_assert ( node_local_id < new_nodes_size );
      

      
      // create a node but pass it the place in memory where it will reside
      // comming from the global data array
      const bool onMesh = true;
      Node* new_node_ptr = new Node ( nodes.getGlobalData(node_local_id), onMesh );
      

      
      // set the node 
      const bool isUpdatable = chset.add_nodes_updatable[last_added_ndix];
      
      new_node_ptr->setParUpdatable( isUpdatable );
      new_node_ptr->setGlobalID( node_global_id );
      

      
      // place the created pointer in the node DataHandle
      nodes[node_local_id] = new_node_ptr;
      
      // put the coords in the added node
      // verify the coords are correctly sized
      RealVector& new_node_coords = chset.add_nodes_coords[last_added_ndix];
      cf_assert ( new_node_ptr->size() == new_node_coords.size() );
      

      
      *new_node_ptr = new_node_coords;
      

      
      // insert the node in the IndexList
      // this will give him his local ID
      IndexList<Node>::getList().createID(new_node_ptr);      
      
      // check that the local ID provided by the IndexList 
      // also matches our accounting
      cf_assert ( new_node_idx == new_node_ptr->getLocalID() );      
      

      
      // increment the index of added nodes
      // this is just for checking purposes
      ++new_node_idx; 
    }
  }
  
  CFLog ( NOTICE , "NEW Nb Nodes [" << new_nodes_size << "]\n" );
}
  
//////////////////////////////////////////////////////////////////////////////

void ChangeMesh::add_remove_states( MeshChangeSet& chset )
{
  CFAUTOTRACE;
  
  DataHandle < State*, GLOBAL >  states = socket_states.getDataHandle();
  
  
  // some initial sanity checks
  
  const CFuint nb_states_to_remove = chset.remove_states.size();
  cf_assert ( chset.add_states_coords.size()    == chset.add_states_gid.size() );
  cf_assert ( chset.add_states_values.size()    == chset.add_states_gid.size() );
  cf_assert ( chset.add_states_updatable.size() == chset.add_states_gid.size() );
  
  
  /// @TODO note that we are not deleting the states that are removed
  ///       we reuse their memory instead of adding a new one, but those 
  ///       that are left are not deleted
  
  const CFuint nb_states_to_increase = 
  chset.add_states_gid.size() > chset.remove_states.size() ? 
  chset.add_states_gid.size() - chset.remove_states.size() : 0 ;
  
  // remember the old number of states
  const CFuint old_states_size = states.size();
  CFLog ( NOTICE , "OLD Nb States [" << old_states_size << "]\n" );
  
    
  
  // resize the container of states
  states.resize( old_states_size + nb_states_to_increase );
  const CFuint new_states_size = states.size();
  
  
  // correctly size the indexes of the states to add
  chset.added_states_local_idx.resize(chset.add_states_gid.size());
  
    
  
  // loop over the old states to remove states
  // and if needed add new ones in their place reusing their local ID
  
  CFuint last_removed_sidx = 0; // will vary from 0 to nb states to remove
  CFuint last_added_sidx = 0;   // will vary from 0 to nb states to add
  
  for ( CFuint i = 0;  i < old_states_size; ++i )
  {
      
    
    CFLog ( NOTICE, "STATE [" << i << "]\n");
    
    bool state_removed = false;
    
      
    
    // verify that we are not passing ahead of state indexes to remove
    // we assume that remove_states is an ordered vector
    if ( nb_states_to_remove != 0 )
    {    
      cf_assert ( i <= chset.remove_states[last_removed_sidx] );
    }
    
    
    // is this state to be removed?
    if ( last_removed_sidx < nb_states_to_remove && chset.remove_states[last_removed_sidx] == i )
    {
      
  
      
      // state memory is not removed
      state_removed = true;
      
      CFLog ( NOTICE, "Remove state with local index [" << i << "]\n");
      
      // incremet index of states to remove
      ++last_removed_sidx; 
    }
    
      
    
    // reuse the removed index if there are states to add
    if ( state_removed && last_added_sidx < chset.add_states_gid.size() ) 
    {
  
      
      const CFuint state_global_id = chset.add_states_gid[last_added_sidx];
      const CFuint state_local_id  = i;
      
      chset.added_states_local_idx[last_added_sidx] = state_local_id;
      
      
      CFLog ( NOTICE, "Replaced with state with global index [" << state_global_id << "]\n");
      
      const bool isUpdatable = chset.add_states_updatable[last_added_sidx];
      
      states[i]->setGhost( false );
      states[i]->setParUpdatable( isUpdatable );
      states[i]->setGlobalID( state_global_id );
      states[i]->setLocalID ( state_global_id );
      
      // put the values in the added state
      // verify the values are correctly sized
      RealVector& new_state_values = chset.add_states_values[last_added_sidx];
      cf_assert ( states[i]->size() == new_state_values.size() );
      
      *states[i] = new_state_values;
      
      // increment index of added states
      ++last_added_sidx;     
    }
    
      
    
  }
  
    
  
  // check all states where removed
  cf_assert ( last_removed_sidx == chset.remove_states.size() );
  
  // check all states where added
  // if not we need to add more points
  if ( last_added_sidx < chset.add_states_gid.size() )
  {
    
    
    
    // start from the last state index
    CFuint new_state_idx = old_states_size;
    
    // loop over the remaining states to add
    // and grow the states memory
    for ( ; last_added_sidx < chset.add_states_gid.size(); ++last_added_sidx )
    {
      // global id of state to add
      const CFuint state_global_id = chset.add_states_gid[last_added_sidx];
      
      CFLog ( NOTICE, "NEW STATE [" << last_added_sidx << "] with GLOBAL ID [" << state_global_id << "]\n");
      

      
      // get the state_local_id from the global ID
      const CFuint state_local_id  = states.addLocalPoint ( state_global_id );
      

      
      // check that the returned local ID matches 
      // the new state idx ccording to our accounting
      cf_assert ( new_state_idx == state_local_id );
      
      // we already resize the states to the new size so the index 
      // should be lower that the size
      cf_assert ( state_local_id < new_states_size );
      

      
      // create a state but pass it the place in memory where it will reside
      // comming from the global data array
      State* new_state_ptr = new State ( states.getGlobalData(state_local_id) );
      

      
      // set the state 
      const bool isUpdatable = chset.add_states_updatable[last_added_sidx];
      
      new_state_ptr->setGhost( false );
      new_state_ptr->setParUpdatable( isUpdatable );
      new_state_ptr->setGlobalID( state_global_id );
      

      
      // place the created pointer in the state DataHandle
      states[state_local_id] = new_state_ptr;
      
      // put the values in the added state
      // verify the values are correctly sized
      RealVector& new_state_values = chset.add_states_values[last_added_sidx];
      cf_assert ( new_state_ptr->size() == new_state_values.size() );
      
      
      /// add the states coordinates ( false -> not owned by the mesh, owned by the state ) 
      Node* new_state_node = new Node ( chset.add_states_coords[last_added_sidx], false );
      
      new_state_ptr->setSpaceCoordinates ( new_state_node );
      
      

      
      *new_state_ptr = new_state_values;
      

      
      // insert the state in the IndexList
      // this will give him his local ID
      IndexList<State>::getList().createID(new_state_ptr);      
      
      // check that the local ID provided by the IndexList 
      // also matches our accounting
      cf_assert ( new_state_idx == new_state_ptr->getLocalID() );      
      

      
      // increment the index of added states
      // this is just for checking purposes
      ++new_state_idx; 
    }
  }
  
  CFLog ( NOTICE , "NEW Nb States [" << new_states_size << "]\n" );
}  

//////////////////////////////////////////////////////////////////////////////
    
void ChangeMesh::add_remove_geoents( MeshChangeSet& chset )     
{    
  SafePtr < std::vector<ElementTypeData > > element_types = MeshDataStack::getActive()->getElementTypeData();
  
  const CFuint nb_elem_types = element_types->size();
  
  for (CFuint iType = 0; iType < nb_elem_types; ++iType)
  {
    // const CFuint nb_states_type = (*element_types)[iType].getNbStates();

    CFLog ( NOTICE, "++++++ Element Type [" << iType << "]\n" );

  }


#if 0
  


  // mesh data builder creates data that are IO-dependent
  meshDataBuilder->computeGeoTypeInfo();
  meshDataBuilder->createTopologicalRegionSets();
  
  // set some global mesh values useful in many places
  // meshDataBuilder->setMaxGlobalInfo();
  meshDataBuilder->setMaxNbStatesInCell();
  meshDataBuilder->setMaxNbNodesInCell();
  meshDataBuilder->setMaxNbFacesInCell();
  
  // meshDataBuilder->releaseMemory();
#endif

}

//////////////////////////////////////////////////////////////////////////////

void ChangeMesh::execute()
{
  CFAUTOTRACE;
	
	const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();

//--------------------------------------------------------------------------------------------

  // DEFINE THE CHANGE TO THE MESH
  
	MeshChangeSet chset;
	
  {
    
    // REMOVE NODES
    chset.remove_nodes.resize(2);
    chset.remove_nodes[0] = 2;
    chset.remove_nodes[1] = 6;
    
    // ADD NODES
    
    chset.add_nodes_gid.resize(6);
    chset.add_nodes_gid[0] = 9;
    chset.add_nodes_gid[1] = 10;
    chset.add_nodes_gid[2] = 11;
    chset.add_nodes_gid[3] = 12;
    chset.add_nodes_gid[4] = 13;
    chset.add_nodes_gid[5] = 14;
    
    chset.add_nodes_coords.resize(6);
    for ( CFuint i = 0; i < chset.add_nodes_coords.size(); ++i) chset.add_nodes_coords[i].resize ( dim );
    chset.add_nodes_coords[0][XX] = 1.5; chset.add_nodes_coords[0][YY] = 0.0;
    chset.add_nodes_coords[1][XX] = 2.0; chset.add_nodes_coords[1][YY] = 0.0;
    chset.add_nodes_coords[2][XX] = 1.5; chset.add_nodes_coords[2][YY] = 0.5;
    chset.add_nodes_coords[3][XX] = 2.0; chset.add_nodes_coords[3][YY] = 0.5;
    chset.add_nodes_coords[4][XX] = 1.5; chset.add_nodes_coords[4][YY] = 1.0;
    chset.add_nodes_coords[5][XX] = 2.0; chset.add_nodes_coords[5][YY] = 1.0;

    chset.add_nodes_updatable.resize(6);
    for ( CFuint i = 0; i < chset.add_nodes_updatable.size(); ++i) chset.add_nodes_updatable[i] = true;
    
    
    // REMOVE STATES
    chset.remove_states.resize(2);
    chset.remove_states[0] = 3;
    chset.remove_states[1] = 4;
    
    // ADD STATES

    chset.add_states_gid.resize(6);
    chset.add_states_gid[0] = 9;
    chset.add_states_gid[1] = 10;
    chset.add_states_gid[2] = 11;
    chset.add_states_gid[3] = 12;
    chset.add_states_gid[4] = 13;
    chset.add_states_gid[5] = 14;
    
    chset.add_states_coords.resize(6);
    for ( CFuint i = 0; i < chset.add_states_coords.size(); ++i) chset.add_states_coords[i].resize ( dim );
    chset.add_states_coords[0][XX] = 1.5; chset.add_states_coords[0][YY] = 0.0;
    chset.add_states_coords[1][XX] = 2.0; chset.add_states_coords[1][YY] = 0.0;
    chset.add_states_coords[2][XX] = 1.5; chset.add_states_coords[2][YY] = 0.5;
    chset.add_states_coords[3][XX] = 2.0; chset.add_states_coords[3][YY] = 0.5;
    chset.add_states_coords[4][XX] = 1.5; chset.add_states_coords[4][YY] = 1.0;
    chset.add_states_coords[5][XX] = 2.0; chset.add_states_coords[5][YY] = 1.0;
    
    chset.add_states_values.resize(6);
    for ( CFuint i = 0; i < chset.add_states_values.size(); ++i) chset.add_states_values[i].resize ( nbeq );
    chset.add_states_values[0] = 0.0;
    chset.add_states_values[1] = 0.0;
    chset.add_states_values[2] = 0.0;
    chset.add_states_values[3] = 0.0;
    chset.add_states_values[4] = 0.0;
    chset.add_states_values[5] = 0.0;
    

    chset.add_states_updatable.resize(6);
    for ( CFuint i = 0; i < chset.add_states_updatable.size(); ++i) chset.add_states_updatable[i] = true;
    
  }
                                                                                       
//--------------------------------------------------------------------------------------------
  
	CFLog ( NOTICE , "\n+++ BEGIN ChangeMesh +++++++++++++++++++++++\n" );
  
  add_remove_nodes ( chset );
  
  CFLog ( NOTICE , "\n---------\n" );
  
  add_remove_states ( chset );

  CFLog ( NOTICE , "\n---------\n" );
  
  add_remove_geoents ( chset );
  
	
	CFLog ( NOTICE , "\n+++ END ChangeMesh +++++++++++++++++++++++\n" );
  
//--------------------------------------------------------------------------------------------
	
  Common::SafePtr<EventHandler> event_handler = 
    Environment::CFEnv::getInstance().getEventHandler();

  // this will resetup all existing methods
  // because they might have private data which depends on the new mesh
  Common::Signal::arg_t msg;
  event_handler->call_signal ( "CF_ON_MESH_UPDATE", msg );
}

//////////////////////////////////////////////////////////////////////////////

void ChangeMesh::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
  
  CFLog ( NOTICE, "++++++ ChangeMesh::configure() [" << getName() << "]\n" );

}
    
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ChangeMesh::needsSockets()
{

  CFLog ( NOTICE, "++++++ ChangeMesh::needsSockets() [" << getName() << "]\n" );

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = DataProcessingCom::needsSockets();
  
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
