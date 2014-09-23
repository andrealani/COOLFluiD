// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "Common/ConnectivityTable.hh"
#include "Common/SwapEmpty.hh"
#include "Common/CFLog.hh"
#include "Common/CFMultiMap.hh"

#include "MathTools/RCM.h"

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

/////////////////////////////////////////////////////////////////////////////

//.......... RCM_algorithm begins here.............

void RCM::renumber (ConnectivityTable<CFuint>& cellstate,
		    ConnectivityTable<CFuint>& cellnode,
		    std::valarray <CFuint>& new_id, 
		    const bool useMedianDual)
{
  CFLog(INFO, "RCM::renumber() => useMedianDual [ " << useMedianDual << " ] START\n");
  
  ConnectivityTable<CFuint> nodenode;
  
  if (!useMedianDual) {
    RCM::transformCellNode2NodeNode (cellstate, nodenode);
  }
  else {
    RCM::transformCellNode2NodeNodeMedianDual (cellstate, cellnode, nodenode);
  }
  
  const CFuint nbelems = cellstate.nbRows();
  const CFuint nnodes  = nodenode.nbRows();
  
  if (useMedianDual) {
    cf_always_assert(nbelems == nnodes);
  }
  
  CFuint readcount   = 0; // readcount = "main counter" on new[i]
  CFuint writecount  = 0; // counter of nodes filling newV[]..

  std::valarray < CFuint > newV ( ( CFuint ) 0, nnodes );
  std::valarray < CFuint > newV_R ( ( CFuint ) 0, nnodes ); // Used in Reverse CUTHILL MCKEE
  std::valarray<bool> flag ( false, nnodes );
  
  // ############################################################################
  // #                            MAIN CYCLE                                    #
  // ############################################################################
  
 
  
  CFuint counter=0; //just a counter..
  while ( readcount < nnodes )
  {
    if ( readcount == 0 )
    {
      
      //............. SEARCHING FOR THE STARTING NODE ................
      CFuint least_conn = std::numeric_limits < CFuint >::max();
      CFuint ileast=0;
      
      for ( CFuint i=0; i<nnodes; ++i )
	{
	  const CFuint nbcols = nodenode.nbCols ( i );
	  if ( (  nbcols < least_conn ) && ( !flag[i] ) )
	    {
	      least_conn = nbcols;
	      ileast=i;
	    }
	};
      
      //..insert into NewV[]
      newV[writecount]=ileast;
      flag[ileast]=true;                  //..so far OK
      
      writecount++;
      
      // 			std::cout << "ileast  " <<ileast<< std::endl;
      //..............................................................
      
      //.... inserting "first-level" nodes into newV[] .... (1)
      CFuint mino=0,minn=0;
      CFuint ccount=0;
      
      for ( CFuint i=0; i < nodenode.nbCols ( newV[readcount] ); i++ )
	{
	  mino = std::numeric_limits < CFuint >::max();
	  minn=mino;
	  
	  for ( CFuint j=0; j < nodenode.nbCols ( newV[readcount] ); j++ )
	    {
	      if ( ( nodenode.nbCols ( nodenode ( newV[readcount],j ) ) <= minn ) && ( !flag[nodenode ( newV[readcount],j ) ] ) )
		{
		  newV[writecount+i] = nodenode ( newV[readcount],j );
		  minn = nodenode.nbCols ( nodenode ( newV[readcount],j ) );
		}
	    }
	  if ( minn != mino )
	    {
	      flag[newV[writecount+i]]=true;
	      ccount++;
	    }
	  
	}
      
      writecount= writecount + ccount;
      readcount++;
      //-----------------------------------------------------
      
      
    } // end if ......1ST ITERATION (placing starting node..)
    
    //Now, I have to repeat all the operations done in the first iteration
    // for all the other nodes, but this time, nodes have to be taken from newV[i]....
    
    //.... inserting "first-level" nodes into newV[] .... (main)
    CFuint mino=0,minn=0;
    CFuint ccount=0;
    
    for ( CFuint i=0; i < nodenode.nbCols ( newV[readcount] ); i++ )
      {
	mino = std::numeric_limits < CFuint >::max();
	minn=mino;
	
	for ( CFuint j=0; j < nodenode.nbCols ( newV[readcount] ); j++ )
	  {
	    
	    if ( ( nodenode.nbCols ( nodenode ( newV[readcount],j ) ) <= minn ) && ( !flag[nodenode ( newV[readcount],j ) ] ) )
	      {
		newV[writecount+i] = nodenode ( newV[readcount],j );
		minn = nodenode.nbCols ( nodenode ( newV[readcount],j ) );
	      }
	  }
	if ( minn != mino )
	  {
	    flag[newV[writecount+i]]=true;
	    ccount++;
	  }
      }
    
    writecount= writecount + ccount;
    ++readcount;
    ++counter;
    // 		if ( ( counter % 100 ) == 0 ) { std::cout << counter << "\n"; }
    //...end cycle on newV[i]......}
  }
  
//----------------- REVERSE CUTHILL MCKEE ---------------------
  
for ( CFuint i=0; i<nnodes ; ++i )
    {
      newV_R[nnodes-i-1]=newV[i];
    }
  
  for ( CFuint i=0; i<nnodes; ++i )
    {
      newV[i]=newV_R[i];
    }

//-------------------------------------------------------------
  
  // make sure that new_id has same number of rows as cellstate
  new_id.resize(nnodes);
  // assign new ids
  for ( CFuint i=0;i<nnodes; ++i ) {
    new_id [ newV[i] ] = i;
  }
  
//--------------------------------------------------------------------------
//                      REWRITE THE TABLE nodenode                               -
//--------------------------------------------------------------------------
  for ( CFuint i=0; i<nnodes; ++i )
    {
      const CFuint nbcols = nodenode.nbCols ( i );
      for ( CFuint j=0; j< nbcols; ++j )
	{
	  nodenode ( i,j ) = new_id[nodenode ( i,j ) ];
	}
    }

//--------------------------------------------------------------------------
//                      REWRITE THE input TABLEs                               -
//--------------------------------------------------------------------------
  if (!useMedianDual) {
    // rewrite the cellstate connectivity table
    for ( CFuint i=0; i< nbelems; ++i ) {
      const CFuint nbcols = cellstate.nbCols ( i );
      for ( CFuint j=0; j<nbcols; ++j) {
	cellstate ( i,j ) = new_id[cellstate ( i,j ) ];
      }
    }
  }
  else {
    // cell-state connectivity is left unchanged for simplicity
    // I can still say that cell 0 has state 0, but is defined by different nodes
    // what will have to change is the actual content (solution vector) of the state
    
    // modify cell-node connectivity by switching "rows", namely cell IDs 
    ConnectivityTable<CFuint> cellnodeBkp(cellnode);
    vector<bool> flag(nbelems, false);
    
    for ( CFuint i=0; i< nbelems; ++i ) {
      const CFuint newCellID = new_id[i];
      // move nodes corresponding to former cell i to newCellID
      const CFuint nbcols = cellnode.nbCols ( i );
      for ( CFuint j=0; j<nbcols; ++j) {
	cellnode ( newCellID,j ) = cellnodeBkp( i,j );
      }
      flag[newCellID] = true;
    }
    
    for ( CFuint i=0; i< nbelems; ++i ) {
      if (!flag[i]) CFLog(ERROR, "ERROR: RCM::renumber() => Nodes for cell [" << i << "] have not been updated!\n");
    }
  }
  
  CFLog(INFO, "RCM::renumber() => END\n");
  
  // function ends
}

/////////////////////////////////////////////////////////////////////////////

void RCM::transformCellNode2NodeNode ( const ConnectivityTable<CFuint>& cellnode, ConnectivityTable<CFuint>& nodenode )
{
  CFuint maxN = 0;
  const CFuint nbrows = cellnode.nbRows();
  
  //........counting number of nodes stored in table cellnode .........
  for ( CFuint i = 0; i < nbrows; ++i )
  {
    const CFuint nbcols = cellnode.nbCols ( i );
    for ( CFuint j = 0; j<nbcols; ++j )
    {
      const CFuint nn = cellnode ( i,j );
      if ( nn >= maxN ) { maxN = nn; }
    }
  }

  const CFuint nnodes = maxN+1;
  //.............................................................

  // construct array with nb of neighbours for each node
  // this is over calculated, because nodes shared between elements get accounted more then once..
  std::valarray <CFuint> maxnb_neighb ( nnodes );
  
  for ( CFuint i = 0; i < nbrows; ++i )
  {
    const CFuint nbcols = cellnode.nbCols ( i );
    cf_assert ( nbcols > 0 );
    const CFuint nnei = nbcols;
    for ( CFuint j = 0; j < nbcols; ++j )
    {
      const CFuint nn = cellnode ( i,j );
      maxnb_neighb[nn] = maxnb_neighb[nn] + nnei;
    }
  }

  //.............................................................
  
  //.... data structure to store the true nb of neighbours .......
  std::vector < std::vector<CFuint> > true_neigh ( nnodes );
  for ( CFuint i = 0; i < nnodes; ++i )
  {
    true_neigh[i].reserve ( maxnb_neighb[i] );
  }

  for ( CFuint i = 0; i < nbrows; i++ )
  {
    const CFuint nbcols = cellnode.nbCols ( i );
    for ( CFuint j = 0; j < nbcols; j++ )
    {
      for ( CFuint k = 0; k < nbcols; k++ )
      {
        true_neigh[cellnode ( i,j ) ].push_back ( cellnode ( i,k ) );
      }
    }
  }

  for ( CFuint i = 0; i < nnodes; ++i )
  {
    // sort the vector so we can then remove duplicated ids
    std::sort ( true_neigh[i].begin(), true_neigh[i].end(), std::less<CFuint>() );
    // place duplicated ids in end of vector
    std::vector<CFuint>::iterator last_id = std::unique ( true_neigh[i].begin(),true_neigh[i].end() );
    // remove duplicated id
    true_neigh[i].erase ( last_id,true_neigh[i].end() );
  }
  
  //.............................................................
  
  //...here I create the table nodenode to store State to State connectivity.....
  //......the pattern first......
  
  std::valarray <CFuint> pat_node ( nnodes );
  for ( CFuint n = 0; n < nnodes; ++n )
  {
    pat_node[n] = true_neigh[n].size();
  }
  
  //... and then the table itself..........
  nodenode.resize ( pat_node );
  for ( CFuint n = 0; n < nnodes; n++ )
  {
    const CFuint nneig = nodenode.nbCols ( n );
    for ( CFuint j = 0; j < nneig; j++ )
    {
      nodenode ( n,j ) = true_neigh[n][j];
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

int RCM::read_input ( const std::string& filename, ConnectivityTable<CFuint>& cellnode )
{
  CFuint nb_elems;

  std::ifstream file_in ( filename.c_str() );
  if ( ! file_in.is_open() )
    {
      std::cout << "Could not open file " <<  filename << std::endl;
      throw std::string ("Could not open file " +  filename);
    }
  
  file_in >> nb_elems;
  
  std::valarray <CFuint> pattern ( nb_elems );
  for ( CFuint i = 0; i < nb_elems; ++i ) // read pattern
    {
      file_in >> pattern[i];
    }
  
  cellnode.resize ( pattern );
  file_in >> cellnode;
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

void RCM::print_table (const std::string& filename, 
		       const ConnectivityTable<CFuint>& cellstate, 
		       const ConnectivityTable<CFuint>& cellnode, 
		       const bool useMedianDual)
{
  ConnectivityTable<CFuint> nodenode;
  
  if (!useMedianDual) {
    RCM::transformCellNode2NodeNode (cellstate, nodenode);
  }
  else {
    RCM::transformCellNode2NodeNodeMedianDual (cellstate, cellnode, nodenode);
  }
  
  const CFuint nnodes  = nodenode.nbRows();
  cf_assert(nnodes > 0);
  if (useMedianDual) {cf_assert(nnodes == cellstate.nbRows());}
  
  //-------- write nodes connectivity in "filename"-------------
  ofstream inputTEC ( filename.c_str() );
  if ( !inputTEC.is_open() )
  {
    throw std::string ("Unable to open file '" +  filename + "' to write data");
  }
  
  inputTEC << "x " <<"y " << "\n";
  
  for ( CFuint n = 0; n < nnodes; n++ )
  {
    const CFuint nneig = nodenode.nbCols ( n );
    for ( CFuint j=0; j<nneig; ++j )
    {
      inputTEC <<  n << " " << nodenode ( n,j )  <<"\n";
    }
  }
  inputTEC.close();
  
  CFLog(VERBOSE, "RCM::print_table() => Wrote file " << filename << "\n");
}
  
/////////////////////////////////////////////////////////////////////////////

void RCM::transformCellNode2NodeNodeMedianDual (const ConnectivityTable<CFuint>& cellstate, 
						const ConnectivityTable<CFuint>& cellnode, 
						ConnectivityTable<CFuint>& nodenode)
{
  using namespace std;

  // count number of pairs node-state to be stored in map node2state
  CFuint nbPairsNodeState = 0;
  const CFuint nbCells = cellnode.nbRows();
  for ( CFuint i = 0; i < nbCells; ++i ) {
    nbPairsNodeState += cellnode.nbCols(i);
  }
  
  // allocate the memory for the map node-state
  typedef CFMultiMap<CFuint, CFuint> MapNodeState;
  MapNodeState mapNodeState(nbPairsNodeState);
  
  // at first a multi-map with key = Node* and value = State*
  // referencing Node is created, by looping over all the nodes
  // inside the element-node connectivity
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    const CFuint nbNodesInCell = cellnode.nbCols(iCell);
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      mapNodeState.insert(cellnode(iCell, iNode), iCell);
    }
  }
  mapNodeState.sortKeys();
  
  cf_assert(mapNodeState.getSize() == nbPairsNodeState);
  
  typedef MapNodeState::MapIterator mapIt;
  typedef set<CFuint, less<CFuint> > SetOfNeighbors;
  
  /// number of neighbors per cell
  std::valarray <CFuint> pat_node (nbCells);
  /// temporary stencil storage
  vector<vector<CFuint> > stencil(nbCells);
  
  // Loop is made over all the cells (remember that iCell == iState !!!)
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    cf_assert(iCell == cellstate(iCell,0));
    const CFuint nbNodesInCell = cellnode.nbCols(iCell);
    SetOfNeighbors* neighborList = new SetOfNeighbors();
    // loop over the nodes of the current cell
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      const CFuint nodeID = cellnode(iCell, iNode);
      // Get all the states referencing that node and store all of
      // them except the state inside the current cell (in fact,
      // you are looking for its neighbors)
      bool fo = false;
      pair<mapIt,mapIt> statesRefNode = mapNodeState.find(nodeID, fo);
      cf_assert(fo);
      
      for (mapIt stateInMap = statesRefNode.first;
	   stateInMap != statesRefNode.second;
	   ++stateInMap) {
	const CFuint stateID = stateInMap->second;
	// if the stateID of the neighbor is less than the current one,
	// the pair formed by this state and the neighbor has been
	// already inserted while processing the neighbor
	if (stateID != iCell) {
	  neighborList->insert(SetOfNeighbors::value_type(stateID));
	}
      }
    }
  
    // store the set of neighbor states in neighborStates
    const CFuint nbNeighStates = neighborList->size();
    stencil[iCell].resize(nbNeighStates);
    pat_node[iCell] = nbNeighStates;
    cf_assert(nbNeighStates > 0);
    
    SetOfNeighbors::iterator it;
    CFuint in = 0;
    for (it = neighborList->begin(); it != neighborList->end(); ++it) {
      stencil[iCell][in] = *it;
      in++;
    }
    deletePtr(neighborList);
  }
  
  // fill in the connectivity table state-2-state (cell-2-cell)
  nodenode.resize(pat_node); 
    
  for ( CFuint n = 0; n < nbCells; n++) {
    const CFuint nneig = nodenode.nbCols (n);
    for (CFuint j = 0; j < nneig; j++) {
      nodenode (n,j) = stencil[n][j];
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD
