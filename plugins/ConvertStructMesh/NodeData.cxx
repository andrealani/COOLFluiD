// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NodeData.hh"
#include "Common/SwapEmpty.hh"
#include "Common/CFMultiMap.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////
  
// void NodeData::assignIDs(MapCoord2ID& bcoord2LocalID, 
// 			 const vector<vector<bool> >& isBnode,
// 			 vector<RealMatrix>& xyzInBlock,
// 			 CFuint nbBNodes)
// {    
//   const CFuint nbBlocks = xyzInBlock.size();
//   CFuint totNbNodes = 0;
//   CFuint maxNbNodesInBlock = 0;
//   for (CFuint ib = 0; ib < nbBlocks; ++ib) {
//     const CFuint currNbNodes = xyzInBlock[ib].nbRows();
//     totNbNodes += currNbNodes;
//     maxNbNodesInBlock = max(maxNbNodesInBlock, currNbNodes);
    
//     // resize the storage of the global node IDs
//     _nodeIDs[ib].resize(currNbNodes);
//     // assign a default value of -1
//     _nodeIDs[ib].assign(currNbNodes,-1);
//   }
  
//   typedef MapCoord2ID::MapIterator MapIt;
//   RealMatrix nodeData(totNbNodes, _dim); // a bit overestimated

//   //  CFuint nodeCount = 0;
// //    for (CFuint ib = 0; ib < nbBlocks; ++ib) {
// //      const RealMatrix& xyz = xyzInBlock[ib];
// //      const CFuint nbBlockNodes = xyz.nbRows(); 
// //      for (CFuint i = 0; i < nbBlockNodes; ++i) {
// //        if (_nodeIDs[ib][i] == -1) {
// //  	const CFreal x = xyz(i,XX);
// //  	const CFreal y = xyz(i,YY);
// //  	const CFreal z = xyz(i,ZZ);
// //  	bool isUnique = true;
	
// // 	// look in the nodes that have jb <= ib and j < i  
// //  	for (CFuint jb = 0; jb <= ib && isUnique; ++jb) {
// //  	  const RealMatrix& coord = xyzInBlock[jb];
// //  	  const CFuint jMax = min(i, xyzInBlock[jb].nbRows());
// //  	  for (CFuint j = 0; j < jMax; ++j) {
// //  	    if ((coord(j,XX) == x) && (coord(j,YY) == y) && (coord(j,ZZ) == z)) {
// //  	      isUnique = false;
// //  	      _nodeIDs[ib][i] = _nodeIDs[jb][j];
// //  	      break;
// //  	    }
// //  	  }
// //  	}
	
// //  	if (isUnique) {
// // 	  // cout << "nodeCount = " << nodeCount << endl;
// //  	  for (CFuint iDim = 0; iDim < _dim; ++iDim) {
// //  	    nodeData(nodeCount, iDim) = xyz(i,iDim); 
// //  	    _nodeIDs[ib][i] = nodeCount;
// //  	  }
// //  	  nodeCount++;
// //  	}
// //        }
// //      }
// //    }
  
  
//   CFuint nodeCount = 0;
//   for (CFuint ib = 0; ib < nbBlocks; ++ib) {
//     const RealMatrix& xyz = xyzInBlock[ib];
//     const CFuint nbBlockNodes = xyz.nbRows(); 
//     for (CFuint i = 0; i < nbBlockNodes; ++i) {
//       if (_nodeIDs[ib][i] == -1) {
// 	if (!isBnode[ib][i]) {
// 	  for (CFuint iDim =0; iDim < _dim; ++iDim) {
// 	    nodeData(nodeCount, iDim) = xyz(i,iDim); 
// 	  }
// 	  _nodeIDs[ib][i] = nodeCount++;
// 	}
// 	else {
// 	  // if the boundary node hasn't got an ID assigned yet,
// 	  //	assign it
// 	  // if (_nodeIDs[ib][i] == -1) {
// 	  bool found = false;
// 	  pair<MapIt, MapIt> bnodes = bcoord2LocalID.
// 	    find(xyz(i,XX), xyz(i,YY), xyz(i,ZZ), found);
	  
// 	  if (!found) {
// 	    cout << "NOT FOUND node (x,y,z) = " << 
// 	      xyz(i,XX) << " " << xyz(i,YY) << " " << xyz(i,ZZ) << endl;
// 	    abort();
// 	  }
	  
// 	  for (CFuint iDim =0; iDim < _dim; ++iDim) {
// 	    nodeData(nodeCount, iDim) = xyz(i,iDim); 
// 	  }
	  
// 	  for (MapIt it = bnodes.first; it != bnodes.second; ++it) {
// 	    CFuint* ptr = it->fourth;
// 	    const CFuint blockID   = ptr[0];
// 	    const CFuint locNodeID = ptr[1];
// 	    _nodeIDs[blockID][locNodeID] = nodeCount;
// 	  }
// 	  nodeCount++;
// 	}
//       }
//     }
//   }
  
//   cf_assert(nodeCount < totNbNodes);

//   SwapEmpty(xyzInBlock);
//   bcoord2LocalID.clear();
  
//   // store all the nodes coordinates without overhead in memory 
//   _allNodes.resize(nodeCount, _dim);
//   for (CFuint i = 0; i < _allNodes.size(); ++i) {
//     _allNodes[i] = nodeData[i];
//   }
  
//   cout << "totNbNodes = " << totNbNodes << endl;
//   cout << "nb unique nodes = " << nodeCount << endl;
// }

void NodeData::assignIDs(MapCoord2ID& bcoord2LocalID, 
			 const vector<vector<bool> >& isBnode,
			 vector<RealMatrix>& xyzInBlock,
			 CFuint nbBNodes)
{    
  const CFuint nbBlocks = xyzInBlock.size();
  CFuint totNbNodes = 0;
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    const CFuint currNbNodes = xyzInBlock[ib].nbRows();
    totNbNodes += currNbNodes;
    
    // resize the storage of the global node IDs
    _nodeIDs[ib].resize(currNbNodes);
    // assign a default value of -1
    _nodeIDs[ib].assign(currNbNodes,-1);
  }
    
  CFMultiMap<Node3D, int*> mapNode3DToID(totNbNodes);
  typedef CFMultiMap<Node3D, int*>::MapIterator MapIt;
  
  Node3D node(CFNULL);
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    const CFuint currNbNodes = xyzInBlock[ib].nbRows();
    for(CFuint iNode = 0; iNode < currNbNodes; ++iNode) {
      node.reset(&xyzInBlock[ib](iNode,XX));
      mapNode3DToID.insert(node, &_nodeIDs[ib][iNode]);
    }
  }
  mapNode3DToID.sortKeys();
  
  RealMatrix nodeData(totNbNodes, _dim); // a bit overestimated
  
  CFuint nodeCount = 0;
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    RealMatrix& xyz = xyzInBlock[ib];
    const CFuint nbBlockNodes = xyz.nbRows(); 
    for (CFuint i = 0; i < nbBlockNodes; ++i) {
      //    for (CFuint iDim = 0; iDim < _dim; ++iDim) {
      // 	nodeData(nodeCount, iDim) = xyz(i,iDim); 
      //       }
      //       _nodeIDs[ib][i] = nodeCount;
      //       nodeCount++;
      
      // if the corresponding node doesn't have an ID yet, 
      // it will get one and all the other nodes equal to it 
      // will get the same value
      if (_nodeIDs[ib][i] == -1) {
	node.reset(&xyz(i,XX));
	
	bool fo = false;
	pair<MapIt, MapIt> nodeList = mapNode3DToID.find(node, fo);
	cf_assert(fo);
	
	//	cout << "matching nodes = " << endl;
	for (MapIt it = nodeList.first; it != nodeList.second; ++it) {
	  if (*it->second != -1) { 
	    cout << "node ID should be -1" << endl;
	    abort();
	  }
	  *it->second = nodeCount;
	  // cout.precision(16);
	  // cout << it->first;
	}
	//cout << endl;
	
	for (CFuint iDim = 0; iDim < _dim; ++iDim) {
	  nodeData(nodeCount, iDim) = xyz(i,iDim); 
	}
	nodeCount++;
      }
      
    }
  }
  
  cf_assert(nodeCount < totNbNodes);
  
  SwapEmpty(xyzInBlock);
  bcoord2LocalID.clear();
  
  // store all the nodes coordinates without overhead in memory 
  _allNodes.resize(nodeCount, _dim);
  for (CFuint i = 0; i < _allNodes.size(); ++i) {
    _allNodes[i] = nodeData[i];
  } 
  
  cout << "totNbNodes = " << totNbNodes << endl;
  cout << "nb unique nodes = " << nodeCount << endl;
}

//////////////////////////////////////////////////////////////////////////////

void NodeData::assignIDs(vector<Block>& blocks, 
			 const vector<vector<bool> >& isBnode,
			 vector<RealMatrix>& xyzInBlock)
{    
  cout << "NodeData::assignIDs() => start" << endl;

  // initialize the _nodeIDs storage
  const CFuint nbBlocks = xyzInBlock.size();
  CFuint totNbNodes = 0;
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    const CFuint currNbNodes = xyzInBlock[ib].nbRows();
    totNbNodes += currNbNodes;
    
    // resize the storage of the global node IDs
    _nodeIDs[ib].resize(currNbNodes);
    // assign a default value of -1
    _nodeIDs[ib].assign(currNbNodes,-1); 
  }
  
  RealMatrix nodeData(totNbNodes, _dim); // a bit overestimated
  CFuint nodeCount = 0;
  
  // assign an ID only to the internal nodes
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    const RealMatrix& xyz = xyzInBlock[ib];
    const CFuint currNbNodes = xyzInBlock[ib].nbRows();
    for (CFuint i = 0; i < currNbNodes; ++i) {
      if (!isBnode[ib][i]) {
	_nodeIDs[ib][i] = nodeCount;
	
	for (CFuint iDim = 0; iDim < _dim; ++iDim) {
	  nodeData(nodeCount, iDim) = xyz(i,iDim); 
	}
	nodeCount++;
      }
    }
  }
  
  cout << "internal volume nodes = " <<  nodeCount << endl;
  cout << "total nb nodes        = " <<  totNbNodes << endl;
  
  const CFuint nbBlockFaces = 6;
  CFuint sizeEdgeNodes = 0;
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    for (CFuint iFace = 0; iFace < nbBlockFaces; ++iFace) {
      sizeEdgeNodes += blocks[ib].face[iFace].edgeNodes.size();
    }
  }
  if (sizeEdgeNodes == 0) {
    cout << "ERROR: sizeEdgeNodes == 0!!"<< endl;
    abort();
  }
  
  cout << "max size edge nodes = " << sizeEdgeNodes << endl;
  CFMultiMap<Node3D, int*> mapNode3DToID(sizeEdgeNodes);
  typedef CFMultiMap<Node3D, int*>::MapIterator MapIt;
  Node3D node(CFNULL);
  Node3D tmpNode(CFNULL);
  CFuint countEdgeNodes = 0;
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    RealMatrix& xyz = xyzInBlock[ib];
    for (CFuint iFace = 0; iFace < nbBlockFaces; ++iFace) {
      // if the current face in the block is not flagged
      // assign IDs to the corresponding face internal nodes
      // within the face 
      const vector<int>& internalNodes = blocks[ib].face[iFace].internalNodes;
      const CFuint nbIntNodes = internalNodes.size();
      
      CFMap<Node3D, CFuint> mapIntNode3DToID;
      if (!blocks[ib].face[iFace].skipThisFace) {
	vector<CFuint> vID(nbIntNodes);
	
	for (CFuint in = 0; in < nbIntNodes; ++in) {
	  const CFuint currNode = internalNodes[in];
	  _nodeIDs[ib][currNode] = nodeCount;
	  
	  if (isBnode[ib][currNode] == false) {
	    cout << "ERROR: isBnode[ib][currNode] == false!" << endl;
	    abort();
	  }
	  
	  for (CFuint iDim = 0; iDim < _dim; ++iDim) {
	    nodeData(nodeCount, iDim) = xyz(currNode,iDim); 
	  }
	  
	  node.reset(&xyz(currNode,XX));
	  mapIntNode3DToID.insert(node, nodeCount);
	  
	  vID[in] = nodeCount;
	  
	  nodeCount++;
	}
	mapIntNode3DToID.sortKeys();
	
	// if the block face is not an external boundary face,
	// flag the corresponding face in the neighbor block
	// to be ignored: this will ensure the assignment of a unique ID 
	// for the internal face nodes
	if (blocks[ib].face[iFace].bcID >= 0) {
	  const CFuint neighBlockID = blocks[ib].face[iFace].neighBlock - 1;
	  const CFuint neighFaceID  = blocks[ib].face[iFace].neighFace - 1;
	  blocks[neighBlockID].face[neighFaceID].skipThisFace = true;
	  
	  const vector<int>& neighborNodes = blocks[neighBlockID].face[neighFaceID].internalNodes;
	  const CFuint nbNeighNodes = neighborNodes.size();
	  if (nbNeighNodes != nbIntNodes) {
	    cout << "nbNeighNodes != nbIntNodes" << endl;
	    abort();
	  }
	  
	  // set the same node ID for the face internal nodes of the 
	  // corresponding face in the neighbor block
	  RealMatrix& neighxyz = xyzInBlock[neighBlockID];
	  for (CFuint in = 0; in < nbNeighNodes; ++in) {
	    const CFuint currNode = neighborNodes[in];
	    node.reset(&neighxyz(currNode,XX));
	    
	    bool found = false;
	    for (CFuint ii = 0; ii < nbIntNodes; ++ii) {
	      const CFuint currN = internalNodes[ii];
	      tmpNode.reset(&xyz(currN,XX));
	     
	      if (!(tmpNode != node)) {
		_nodeIDs[neighBlockID][currNode] = vID[ii];
		found = true;
		break;
	      }
	    }
	    if (!found) { 
	      cout << "### blockID = " << ib << ", face = " << iFace << " doesn't match with ";
	      cout << "neighbor blockID = " << neighBlockID << ", face = " << neighFaceID << endl;
	      abort();
	      //    // -------------- remove this ----------------------- //
	      // 	      ofstream* fout1 = new ofstream("face1.plt");  
	      	      
	      // reorder neighbor internal nodes according to the i,j,k convention 
	      // for the current internal nodes  
	      //  if (((iFace == 0 || iFace == 5 || iFace == 2 || iFace == 3) && 
	      // 		   (neighFaceID == 1 || neighFaceID == 4)) ||
	      // 		  ((neighFaceID == 0 || neighFaceID == 5 || neighFaceID == 2 || neighFaceID == 3) && 
	      // 		   (iFace == 1 || iFace == 4))) {
	      // 		vector<CFuint> neworder(nbNeighNodes);
	      // 		CFuint n1 = 0;
	      // 		CFuint n2 = 0;
	      // 		switch(neighFaceID) {
	      // 		case(0):
	      // 		  n1 = _nx[neighBlockID];
	      // 		  n2 = _ny[neighBlockID];
	      // 		  break;
	      // 		case(1):
	      // 		  n1 = _nx[neighBlockID];
	      // 		  n2 = _nz[neighBlockID];
	      // 		  break;
	      // 		case(2):
	      // 		  n1 = _nz[neighBlockID];
	      // 		  n2 = _ny[neighBlockID];
	      // 		  break;
	      // 		case(3):
	      // 		  n1 = _nz[neighBlockID];
	      // 		  n2 = _ny[neighBlockID];
	      // 		  break;
	      // 		case(4):
	      // 		  n1 = _nx[neighBlockID];
	      // 		  n2 = _nz[neighBlockID];
	      // 		  break;
	      // 		case(5):
	      // 		  n1 = _nx[neighBlockID];
	      // 		  n2 = _ny[neighBlockID];
	      // 		  break;
	      // 		}
	      
	      // 	      }
	      // 	      else {
	      
		
	      
	      // 	      }
	      
	      // 	      break;
	      // 	      cout << "tmpNodes:" << endl;
	      // 	      for (CFuint ii = 0; ii < nbIntNodes; ++ii) {
	      // 		const CFuint currN = internalNodes[ii];
	      // 		tmpNode.reset(&xyz(currN,XX));
	      // 		*fout1 << tmpNode;
	      // 	      }
	      
	      // 	      delete fout1;
	      // 	      ofstream* fout2 = new ofstream("face2.plt");  
	      
	      // 	      cout << "neighborNodes:" << endl;
	      // 	      for (CFuint in = 0; in < nbNeighNodes; ++in) {
	      // 		const CFuint currNode = neighborNodes[in];
	      // 		tmpNode.reset(&neighxyz(currNode,XX));
	      // 		*fout2 << tmpNode;
	      // 	      }
	      // 	      delete fout2;
	      
	      //cout << "###target NODE =" << node << " not found!" << endl; abort();
	    }

	    
// 	    try {
// 	      _nodeIDs[neighBlockID][currNode] = mapIntNode3DToID.find(node);
// 	    }
// 	    catch (NoSuchValueException& e) {
// 	      cout << "UNMATCHED = " << node << endl;
	      
// 	      for (CFuint n = 0; n < nbIntNodes; ++n) {
// 		const CFuint currN = internalNodes[n];
// 		tmpNode.reset(&xyz(currN,XX));
		
// 		cout << "node to match = " <<  setw(16) << fixed << setprecision(8) <<  tmpNode;		
// 		cout << "diffX = " <<  setw(16) << fixed << setprecision(8)
// 		     <<  fabs(tmpNode[XX]-node[XX]) << " " 
// 		     <<  fabs(tmpNode[YY]-node[YY]) << " " 
// 		     <<  fabs(tmpNode[ZZ]-node[ZZ]) <<  endl << endl;
// 	      }
// 	      abort();	      
// 	    }
	    
	  }
	}
      }
      
      const vector<int>& edgeNodes = blocks[ib].face[iFace].edgeNodes;
      const CFuint nbEdgeNodes = edgeNodes.size();
      for (CFuint in = 0; in < nbEdgeNodes; ++in) {
	const CFuint currNode = edgeNodes[in];  
	
	if (isBnode[ib][currNode] == false) {
	  cout << "ERROR: isBnode[ib][currNode] == false!" << endl;
	  abort();
	}
	
	node.reset(&xyz(currNode,XX));
	mapNode3DToID.insert(node, &_nodeIDs[ib][currNode]);
	countEdgeNodes++;
      }
    }
  }
  mapNode3DToID.sortKeys();
  
  cout << "nodeCount after map = " << countEdgeNodes  << endl;
  cout << "countEdgeNodes after map = " << countEdgeNodes  << endl;
  
  CFuint countNotAssigned = 0;
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    for (CFuint i = 0; i < _nodeIDs[ib].size(); ++i) {
      if (_nodeIDs[ib][i] == -1) {
	countNotAssigned++;
      }
    }
  }
  cout << "countNotAssigned = " << countNotAssigned  << endl;
  cout << "nodeCount = " << nodeCount  << endl;
  
  countEdgeNodes = 0;
  for (CFuint ib = 0; ib < nbBlocks; ++ib) {
    RealMatrix& xyz = xyzInBlock[ib];
    const CFuint nbBlockNodes = xyz.nbRows(); 
    for (CFuint i = 0; i < nbBlockNodes; ++i) {
      // if the corresponding node doesn't have an ID yet, 
      // it will get one and all the other nodes equal to it 
      // will get the same value
      if (_nodeIDs[ib][i] == -1) {
	node.reset(&xyz(i,XX));

	bool fo = false;
	pair<MapIt, MapIt> nodeList = mapNode3DToID.find(node, fo);
	cf_assert(fo);
	
	bool anyNodeFound = false;
	for (MapIt it = nodeList.first; it != nodeList.second; ++it) {
	  if (*it->second != -1) {
	    cout << "isBnode " << isBnode[ib][i] << endl;
	    cout << "node ID should be -1" << endl;
	    abort();
	  }
	  *it->second = nodeCount;
	  countEdgeNodes++;
	  anyNodeFound = true;
	  
	  //  if (*it->second == -1) {
	  // 	    *it->second = nodeCount;
	  // 	    anyNodeFound = true;
	  // 	  }
	  // 	  countEdgeNodes++;
	}
	if (anyNodeFound) {
	  for (CFuint iDim = 0; iDim < _dim; ++iDim) {
	    nodeData(nodeCount, iDim) = xyz(i,iDim); 
	  }
	  nodeCount++;
	}
	else {
	  cout << "ERROR: no matching node found"<< endl;
	  abort();
	}
      }
    }
  }  
  
  if (countEdgeNodes > sizeEdgeNodes) {
    cout << "ERROR: countEdgeNodes > sizeEdgeNodes!" << endl;
    abort();
  }
  cf_assert(nodeCount < totNbNodes);
  
  SwapEmpty(xyzInBlock);
  
  // store all the nodes coordinates without overhead in memory 
  _allNodes.resize(nodeCount, _dim);
  for (CFuint i = 0; i < _allNodes.size(); ++i) {
    _allNodes[i] = nodeData[i];
  } 
  
  cout << "totNbNodes = " << totNbNodes << endl;
  cout << "nb unique nodes = " << nodeCount << endl;
  
  cout << "NodeData::assignIDs() => end" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void NodeData::assignDefaultIDs(CFuint iBlock, CFuint nbNodes)
{
  _nodeIDs[iBlock].resize(nbNodes);
  for (CFuint i = 0; i < nbNodes; ++i) {
    _nodeIDs[iBlock][i] = i;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
