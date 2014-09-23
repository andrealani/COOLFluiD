// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <iostream>
#include <iomanip>

#include "ConvertStructMesh/ConvertGridProMesh.hh"
#include "MathTools/CFMat.hh"
#include "MathTools/CFVec.hh"
#include "Common/CFMap.hh"
#include "ConvertStructMesh/ConvertStructMeshApp.hh"
#include "Environment/ObjectProvider.hh"

//#define BLOCKID 0
//#define CHANGE_PRECISION 0

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ConvertGridProMesh, ConvertStructMesh,
			    ConvertStructMeshAppModule, 1>
convertGridProMeshProv("GridPro");

//////////////////////////////////////////////////////////////////////////////

ConvertGridProMesh::ConvertGridProMesh(const std::string& fileName) :
  ConvertStructMesh(fileName)
{
  _elemIDs.resize(4);
  FaceData::setSize(DIM_3D);
}

//////////////////////////////////////////////////////////////////////////////

ConvertGridProMesh::~ConvertGridProMesh()
{
}

//////////////////////////////////////////////////////////////////////////////

// void ConvertGridProMesh::readAndBuildNodes()
// {
//   cout << "ConvertGridProMesh::readAndBuildNodes() => start"<< endl;

//  ifstream fin(getGlobalMeshFile().c_str());

//   fin >> _nbBlocks;
//   cout << "nbBlocks = " << _nbBlocks << endl;

//   /// nodes per block
//   vector<RealMatrix> xyzInBlock(_nbBlocks);
//   _nodes.setNbBlocks(_nbBlocks);
//   _nodes.setDim(DIM_3D);

//   vector<CFuint>& nxvec = _nodes.getNxVec();
//   vector<CFuint>& nyvec = _nodes.getNyVec();
//   vector<CFuint>& nzvec = _nodes.getNzVec();
//   nxvec.resize(_nbBlocks);
//   nyvec.resize(_nbBlocks);
//   nzvec.resize(_nbBlocks);

//   vector<vector<bool> > isBnode(_nbBlocks);

// #ifdef BLOCKID
//   ofstream fout("blocks.vol");
//   fout << 1 << endl;
// #endif

// #ifdef CHANGE_PRECISION
//   ofstream fout("expert10.vol");
//   fout << _nbBlocks << endl;
// #endif

//   CFuint nbBNodes = 0;
//   _nbElements = 0;
//   for (CFuint ib = 0; ib < _nbBlocks; ++ib) {
//     std::string tmp; // /vol3d
//     fin >> tmp;

// #ifdef BLOCKID
//     if (ib == BLOCKID) {
//       fout << tmp << endl;
//     }
// #endif

// #ifdef CHANGE_PRECISION
//     fout << tmp << endl;
// #endif

//     // maximum number of nodes in x,y,z directions
//     fin >> nxvec[ib] >> nyvec[ib] >> nzvec[ib];

// #ifdef BLOCKID
//     if (ib == BLOCKID) {
//       cout << nxvec[ib] << " " << nyvec[ib] << " " << nzvec[ib] << endl;
//       fout << nxvec[ib] << " " << nyvec[ib] << " " << nzvec[ib] << endl;
//     }
// #endif

// #ifdef CHANGE_PRECISION
//     fout << nxvec[ib] << " " << nyvec[ib] << " " << nzvec[ib] << endl;
// #endif

//     const CFuint totNbNodes = nxvec[ib]*nyvec[ib]*nzvec[ib];
//     xyzInBlock[ib].resize(totNbNodes, DIM_3D);
//     isBnode[ib].resize(totNbNodes, false);

//     const CFuint nx = nxvec[ib];
//     const CFuint ny = nyvec[ib];
//     const CFuint nz = nzvec[ib];
//     const CFuint imax = nx-1;
//     const CFuint jmax = ny-1;
//     const CFuint kmax = nz-1;
//     _nbElements += imax*jmax*kmax;

//     CFuint iNode = 0;
//     for (CFuint k = 0; k < nz; ++k) {
//       for (CFuint j = 0; j < ny; ++j) {
//         for (CFuint i = 0; i < nx; ++i, ++iNode) {
//           cf_assert(iNode < totNbNodes);
// 	  fin >> xyzInBlock[ib](iNode,XX)
// 	      >> xyzInBlock[ib](iNode,YY)
// 	      >> xyzInBlock[ib](iNode,ZZ);

// #ifdef BLOCKID
// 	  if (ib == BLOCKID) {
// 	    fout << setw(16) << fixed << setprecision(12)
// 		 << xyzInBlock[ib](iNode,XX) << " "
// 		 << xyzInBlock[ib](iNode,YY) << " "
// 		 << xyzInBlock[ib](iNode,ZZ) << endl;
// 	  }
// #endif

// #ifdef CHANGE_PRECISION
// 	//   long int tmpX = static_cast<long int>(xyzInBlock[ib](iNode,XX)*10000000);
// // 	  long int tmpY = static_cast<long int>(xyzInBlock[ib](iNode,YY)*10000000);
// // 	  long int tmpZ = static_cast<long int>(xyzInBlock[ib](iNode,ZZ)*10000000);
// // 	  const CFreal x = tmpX/10000000.;
// // 	  const CFreal y = tmpY/10000000.;
// // 	  const CFreal z = tmpZ/10000000.;

// // 	  fout << setw(16) << fixed << setprecision(12)
// // 	       << x << " " << y << " " << z << endl;

// 	  fout << setw(16) << fixed << setprecision(12)
// 	       << xyzInBlock[ib](iNode,XX) << " "
// 	       << xyzInBlock[ib](iNode,YY) << " "
// 	       << xyzInBlock[ib](iNode,ZZ) << endl;
// #endif

// 	  // store separately the boundary nodes also in a map
// 	  if ((i == 0) || (i == imax) ||
// 	      (j == 0) || (j == jmax) ||
// 	      (k == 0) || (k == kmax)) {
// 	    isBnode[ib][iNode] = true;
// 	    nbBNodes++;
// 	  }
// 	}
//       }
//     }
//   }

// #ifdef CHANGE_PRECISION
//   cout << "file with changed precision has been written"<< endl;
//   abort();
// #endif

//   // now we can allocate the map3D with the right size
//   MapCoord2ID bcoord2LocalID(nbBNodes);
//   // array storing the block ID and the local node ID in the block
//   // of a boundary node
//   vector<CFuint> bNodeBlockLocID;
//   bNodeBlockLocID.reserve(nbBNodes*2);

//   CFuint counter = 0;
//   for (CFuint ib = 0; ib < _nbBlocks; ++ib) {
//     const CFuint totNbNodes = nxvec[ib]*nyvec[ib]*nzvec[ib];
//     const CFuint nx = nxvec[ib];
//     const CFuint ny = nyvec[ib];
//     const CFuint nz = nzvec[ib];
//     const CFuint imax = nx-1;
//     const CFuint jmax = ny-1;
//     const CFuint kmax = nz-1;

//     CFuint iNode = 0;
//     for (CFuint k = 0; k < nz; ++k) {
//       for (CFuint j = 0; j < ny; ++j) {
//         for (CFuint i = 0; i < nx; ++i, ++iNode) {
//           cf_assert(iNode < totNbNodes);

// 	  // store separately the boundary nodes also in a map
// 	  if ((i == 0) || (i == imax) ||
// 	      (j == 0) || (j == jmax) ||
// 	      (k == 0) || (k == kmax)) {
// 	    cf_assert(isBnode[ib][iNode]);

// 	    bNodeBlockLocID.push_back(ib);
// 	    bNodeBlockLocID.push_back(iNode);
// 	    // store a mapping (x,y,z) -> pointer to the block
// 	    // and local node ID in the block
// 	    cf_assert(counter < nbBNodes);

// 	    bcoord2LocalID.insert(xyzInBlock[ib](iNode,XX),
// 				  xyzInBlock[ib](iNode,YY),
// 				  xyzInBlock[ib](iNode,ZZ),
// 				  &bNodeBlockLocID[counter*2]);
// 	    counter++;
// 	  }
// 	}
//       }
//     }
//   }

//   cf_assert(counter == nbBNodes);

//   bcoord2LocalID.sortKeys();

//   _nodes.assignIDs(bcoord2LocalID, isBnode, xyzInBlock, counter);

//   cout << "ConvertGridProMesh::readAndBuildNodes() => end"<< endl;
// }

void ConvertGridProMesh::readAndBuildNodes()
{
 cout << "ConvertGridProMesh::readAndBuildNodes() => start"<< endl;

 ifstream fin(getGlobalMeshFile().c_str());

 fin >> _nbBlocks;
 cout << "nbBlocks = " << _nbBlocks << endl;

 // nodes per block
 vector<RealMatrix> xyzInBlock(_nbBlocks);
 _nodes.setNbBlocks(_nbBlocks);
 _nodes.setDim(DIM_3D);

 vector<CFuint>& nxvec = _nodes.getNxVec();
 vector<CFuint>& nyvec = _nodes.getNyVec();
 vector<CFuint>& nzvec = _nodes.getNzVec();
 nxvec.resize(_nbBlocks);
 nyvec.resize(_nbBlocks);
 nzvec.resize(_nbBlocks);

 vector<vector<bool> > isBnode(_nbBlocks);

#ifdef BLOCKID
  ofstream fout("blocks.vol");
  fout << 1 << endl;
#endif

#ifdef CHANGE_PRECISION
  ofstream fout("expert5.vol");
  fout << _nbBlocks << endl;
#endif

  CFuint nbBNodes = 0;
  CFuint nbBInternalNodes = 0;
  CFuint nbBEdgeNodes = 0;

  _nbElements = 0;
  for (CFuint ib = 0; ib < _nbBlocks; ++ib) {
    std::string tmp; // /vol3d
    fin >> tmp;

#ifdef BLOCKID
    if (ib == BLOCKID) {
      fout << tmp << endl;
    }
#endif

#ifdef CHANGE_PRECISION
    fout << tmp << endl;
#endif

    // maximum number of nodes in x,y,z directions
    fin >> nxvec[ib] >> nyvec[ib] >> nzvec[ib];

#ifdef BLOCKID
    if (ib == BLOCKID) {
      cout << nxvec[ib] << " " << nyvec[ib] << " " << nzvec[ib] << endl;
      fout << nxvec[ib] << " " << nyvec[ib] << " " << nzvec[ib] << endl;
    }
#endif

#ifdef CHANGE_PRECISION
    fout << nxvec[ib] << " " << nyvec[ib] << " " << nzvec[ib] << endl;
#endif

    const CFuint totNbNodes = nxvec[ib]*nyvec[ib]*nzvec[ib];
    xyzInBlock[ib].resize(totNbNodes, DIM_3D);
    isBnode[ib].resize(totNbNodes, false);

    const CFuint nx = nxvec[ib];
    const CFuint ny = nyvec[ib];
    const CFuint nz = nzvec[ib];
    const CFuint imax = nx-1;
    const CFuint jmax = ny-1;
    const CFuint kmax = nz-1;
    _nbElements += imax*jmax*kmax;

    CFuint iNode = 0;
    for (CFuint k = 0; k < nz; ++k) {
      for (CFuint j = 0; j < ny; ++j) {
        for (CFuint i = 0; i < nx; ++i, ++iNode) {
          cf_assert(iNode < totNbNodes);
	  fin >> xyzInBlock[ib](iNode,XX)
	      >> xyzInBlock[ib](iNode,YY)
	      >> xyzInBlock[ib](iNode,ZZ);

#ifdef BLOCKID
	  if (ib == BLOCKID) {
	    fout << setw(16) << fixed << setprecision(12)
		 << xyzInBlock[ib](iNode,XX) << " "
		 << xyzInBlock[ib](iNode,YY) << " "
		 << xyzInBlock[ib](iNode,ZZ) << endl;
	  }
#endif

#ifdef CHANGE_PRECISION
	//   long int tmpX = static_cast<long int>(xyzInBlock[ib](iNode,XX)*10000000);
	//   long int tmpY = static_cast<long int>(xyzInBlock[ib](iNode,YY)*10000000);
// 	  long int tmpZ = static_cast<long int>(xyzInBlock[ib](iNode,ZZ)*10000000);
// 	  const CFreal x = tmpX/10000000.;
// 	  const CFreal y = tmpY/10000000.;
// 	  const CFreal z = tmpZ/10000000.;

 	  fout << setw(16) << fixed << setprecision(12)
	       << xyzInBlock[ib](iNode,XX) << " "
	       << xyzInBlock[ib](iNode,YY) << " "
	       << xyzInBlock[ib](iNode,ZZ) << endl;
#endif

	  // store separately the boundary nodes also in a map
	  if ((i == 0) || (i == imax) ||
	      (j == 0) || (j == jmax) ||
	      (k == 0) || (k == kmax)) {
	    isBnode[ib][iNode] = true;
	    nbBNodes++;

	    if (i == 0) {
	      if((j != 0) && (j != jmax) && (k != 0) && (k != kmax)) {
		_blocks[ib].face[2].internalNodes.push_back(iNode);
		nbBInternalNodes++;
	      }
	      else {
		_blocks[ib].face[2].edgeNodes.push_back(iNode);
		nbBEdgeNodes++;
	      }
	      continue;
	    }

	    if (i == imax) {
	      if((j != 0) && (j != jmax) && (k != 0) && (k != kmax)) {
		_blocks[ib].face[3].internalNodes.push_back(iNode);
		nbBInternalNodes++;
	      }
	      else {
		_blocks[ib].face[3].edgeNodes.push_back(iNode);
		nbBEdgeNodes++;
	      }
	      continue;
	    }

	    if (j == 0) {
	      if((i != 0) && (i != imax) && (k != 0) && (k !=kmax)) {
		_blocks[ib].face[1].internalNodes.push_back(iNode);
	     	nbBInternalNodes++;
	      }
	      else {
		_blocks[ib].face[1].edgeNodes.push_back(iNode);
		nbBEdgeNodes++;
	      }
	      continue;
	    }

	    if (j == jmax) {
	      if((i != 0) && (i != imax) && (k != 0) && (k !=kmax)) {
		_blocks[ib].face[4].internalNodes.push_back(iNode);
	    	nbBInternalNodes++;
	      }
	      else {
		_blocks[ib].face[4].edgeNodes.push_back(iNode);
		nbBEdgeNodes++;
	      }
	      continue;
	    }

	    if (k == 0) {
	      if((i != 0) && (i != imax) && (j != 0) && (j !=jmax)) {
		_blocks[ib].face[0].internalNodes.push_back(iNode);
	     	nbBInternalNodes++;
	      }
	      else {
		_blocks[ib].face[0].edgeNodes.push_back(iNode);
		nbBEdgeNodes++;
	      }
	      continue;
	    }

	    if (k == kmax) {
	      if((i != 0) && (i != imax) && (j != 0) && (j !=jmax)) {
		_blocks[ib].face[5].internalNodes.push_back(iNode);
	     	nbBInternalNodes++;
	      }
	      else {
		_blocks[ib].face[5].edgeNodes.push_back(iNode);
		nbBEdgeNodes++;
	      }
	      continue;
	    }
	  }
	}
      }
    }
  }

  cout << "nbBNodes = " << nbBNodes << endl;
  cout << "nbBInternalNodes = " << nbBInternalNodes << endl;
  cout << "nbBEdgeNodes = " << nbBEdgeNodes << endl;

#ifdef CHANGE_PRECISION
  cout << "file with changed precision has been written"<< endl;
  abort();
#endif

  _nodes.assignIDs(_blocks, isBnode, xyzInBlock);

  cout << "ConvertGridProMesh::readAndBuildNodes() => end"<< endl;
}

/////////////////////////////////////////////////////////////

void ConvertGridProMesh::buildElements()
{
  cout << "ConvertGridProMesh::buildElements() => start"<< endl;

  std::string bndFile = _fileName + ".bnd";
  cout << "reading file = " <<  bndFile<< endl;

  vector<CFuint>& nxvec = _nodes.getNxVec();
  vector<CFuint>& nyvec = _nodes.getNyVec();
  vector<CFuint>& nzvec = _nodes.getNzVec();

  ifstream fin(bndFile.c_str());

#ifdef BLOCKID
  ofstream fout("blocks.bnd");
#endif
  typedef CFMultiMap<CFuint,pair<int,int> > BlockFaceMap;
  BlockFaceMap blockToFaceBcID(_nbBlocks*6);

  CFuint blockID = 0;
  pair<int, int> faceBcID;
  CFuint faceCounter = 0;
  bool doRead = true;
  CFuint ib = 0;
  do {
    fin >> blockID >> faceBcID.first >> faceBcID.second;

#ifdef BLOCKID
    if (blockID == BLOCKID+1) {
      fout << 1 << " " << faceBcID.first << " " << faceBcID.second << endl;
      fout << 0 << " " << 0 << " " << 0 << endl;
      fout.close();
    }
#endif

    if ((blockID == 0) && (faceBcID.first == 0) && (faceBcID.second == 0)) {
      doRead = false;
    }
    else {
      const CFuint cfBlockID = blockID - 1;
      blockToFaceBcID.insert(cfBlockID, faceBcID);
      faceCounter += getNbFacesInBlockFace(nxvec[cfBlockID]-1,
					 nyvec[cfBlockID]-1,
					   nzvec[cfBlockID]-1,
					   faceBcID.first);
    }

    ib++;

  } while(doRead);

  fin.close();

#ifdef BLOCKID
  fout.close();
#endif

  blockToFaceBcID.sortKeys();

  cout << bndFile << " read"<< endl;

  _elements.resize(_nbElements,8);
  _bface.reserve(faceCounter);
  _mapBcIDToFace.reserve(faceCounter);

  typedef BlockFaceMap::MapIterator MapIt;

  set<int> bcIDList;

  // build the elements
  CFuint iElem = 0;
  for (CFuint ib = 0; ib < _nbBlocks; ++ib) {
    _imax = nxvec[ib] - 1;
    _jmax = nyvec[ib] - 1;
    _kmax = nzvec[ib] - 1;

    bool blockFound = false;
    pair<MapIt,MapIt> blockFaces = blockToFaceBcID.find(ib, blockFound);

    for (CFuint k = 0; k < _kmax; ++k) {
      for (CFuint j = 0; j < _jmax; ++j) {
        for (CFuint i = 0; i < _imax; ++i, ++iElem) {
	  cf_assert(iElem < _nbElements);

	  _elements(iElem,0) = _nodes.getNodeID(ib,i,j,k);
	  _elements(iElem,1) = _nodes.getNodeID(ib,i+1,j,k);
	  _elements(iElem,2) = _nodes.getNodeID(ib,i+1,j+1,k);
	  _elements(iElem,3) = _nodes.getNodeID(ib,i,j+1,k);
	  _elements(iElem,4) = _nodes.getNodeID(ib,i,j,k+1);
	  _elements(iElem,5) = _nodes.getNodeID(ib,i+1,j,k+1);
	  _elements(iElem,6) = _nodes.getNodeID(ib,i+1,j+1,k+1);
	  _elements(iElem,7) = _nodes.getNodeID(ib,i,j+1,k+1);

	  if (blockFound) {
	    // if the element is a boundary element
	    if ((i == 0 || i == _imax-1) ||
		(j == 0 || j == _jmax-1) ||
		(k == 0 || k == _kmax-1)) {
	      // check if the corresponding face of the element is referenced in the map
	      // if yes, add the face in the FaceData storage
	      for (MapIt f = blockFaces.first; f != blockFaces.second; ++f) {
                const CFuint hexaFaceID = f->second.first;
                const CFuint bcID       = f->second.second;
	        FaceData::Itr bFaceItr = addBoundaryFace(iElem, i, j, k, hexaFaceID, bcID);
		if (!bFaceItr.isNull()) {
	          _mapBcIDToFace.insert(bcID, bFaceItr);
	          bcIDList.insert(bcID);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  _mapBcIDToFace.sortKeys();

  _bcListIDs.reserve(bcIDList.size());
  set<int>::iterator it;
  for (it = bcIDList.begin(); it != bcIDList.end(); ++it) {
    _bcListIDs.push_back(*it);
  }

  if (_bface.getNbFaces() !=  faceCounter) {
    cout << "ERROR : _bface.getNbFaces() != faceCounter" << endl;
    cout << "_bface.getNbFaces() = " << _bface.getNbFaces() << endl;
    cout << "faceCounter = " << faceCounter << endl;
    abort();
  }
  else {
    cout << "nb bfaces detected = " << _bface.getNbFaces() << endl;
    cout << "nb BCs = " << bcIDList.size() << endl;
  }

  cout << _nbElements << " elements built" << endl;
  cout << "ConvertGridProMesh::buildElements() => end"<< endl;
}

//////////////////////////////////////////////////////////////////////////////


CFuint ConvertGridProMesh::getNbFacesInBlockFace(CFuint maxnx, CFuint maxny,
						CFuint maxnz, int hexaFaceID) const
{
  if (hexaFaceID == 1 || hexaFaceID == 6) {
    return maxnx*maxny;
  }

  if (hexaFaceID == 3 || hexaFaceID == 4) {
    return maxny*maxnz;
  }

  if (hexaFaceID == 2 || hexaFaceID == 5) {
    return maxnx*maxnz;
  }

  cout << "ConvertGridProMesh::getNbFacesInBlockFace() => hexaFaceID < 1 || hexaFaceID > 6\n";
  abort();
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

FaceData::Itr ConvertGridProMesh::addBoundaryFace(CFuint iElem, CFuint i, CFuint j,
						 CFuint k, int hexaFaceID, int bcID)
{
  if (k == 0 && hexaFaceID == 1) {
    _elemIDs[0] = _elements(iElem,0);
    _elemIDs[1] = _elements(iElem,3);
    _elemIDs[2] = _elements(iElem,2);
    _elemIDs[3] = _elements(iElem,1);
    return _bface.addData(hexaFaceID,iElem, bcID,_elemIDs);
  }

  if (k == _kmax-1 && hexaFaceID == 6) {
    _elemIDs[0] = _elements(iElem,4);
    _elemIDs[1] = _elements(iElem,5);
    _elemIDs[2] = _elements(iElem,6);
    _elemIDs[3] = _elements(iElem,7);
    return _bface.addData(hexaFaceID,iElem, bcID,_elemIDs);
  }

  if (i == 0 && hexaFaceID == 3) {
    _elemIDs[0] = _elements(iElem,0);
    _elemIDs[1] = _elements(iElem,4);
    _elemIDs[2] = _elements(iElem,7);
    _elemIDs[3] = _elements(iElem,3);
    return _bface.addData(hexaFaceID,iElem, bcID,_elemIDs);
  }

  if (i == _imax-1 && hexaFaceID == 4) {
    _elemIDs[0] = _elements(iElem,1);
    _elemIDs[1] = _elements(iElem,2);
    _elemIDs[2] = _elements(iElem,6);
    _elemIDs[3] = _elements(iElem,5);
    return _bface.addData(hexaFaceID,iElem, bcID,_elemIDs);
  }

  if (j == 0 && hexaFaceID == 2) {
    _elemIDs[0] = _elements(iElem,0);
    _elemIDs[1] = _elements(iElem,1);
    _elemIDs[2] = _elements(iElem,5);
    _elemIDs[3] = _elements(iElem,4);
    return _bface.addData(hexaFaceID,iElem, bcID,_elemIDs);
  }

  if (j == _jmax-1 && hexaFaceID == 5) {
    _elemIDs[0] = _elements(iElem,2);
    _elemIDs[1] = _elements(iElem,3);
    _elemIDs[2] = _elements(iElem,7);
    _elemIDs[3] = _elements(iElem,6);
    return _bface.addData(hexaFaceID,iElem, bcID,_elemIDs);
  }

  return FaceData::Itr(CFNULL);
}

//////////////////////////////////////////////////////////////////////////////

void ConvertGridProMesh::writeBoundaryTecplotFile()
{
  cout << "ConvertGridProMesh::writeBoundaryTecplotFile() => start" << endl;

  typedef CFMultiMap<int,FaceData::Itr>::MapIterator MapIt;

  std::string fileName = _fileName + ".surf.plt";
  ofstream fout(fileName.c_str());

  fout << "TITLE =  Unstructured grid data" << endl;
  fout << "VARIABLES  =  \"x0\" \"x1\" \"x2\" " << endl;

  vector<CFuint> nodeIDsinTRS;
  vector<int>::iterator it;
  CFuint totalNbFaces = 0;

  for (it = _bcListIDs.begin(); it != _bcListIDs.end(); ++it) {
    cout << "Bc ID = " << *it << " ";
    
    bool fo = false;
    pair<MapIt,MapIt> faces = _mapBcIDToFace.find(*it, fo);
    if (!fo) cout << "ConvertGridProMesh::writeBoundaryTecplotFile() => face " << *it << " not found!\n";
    
    CFuint nbFaces = 0;
    for (MapIt f = faces.first; f != faces.second; ++f, ++nbFaces, ++totalNbFaces) {
      int* nodesInFace = f->second.getNodesID();
      for (CFuint i =0; i < 4; ++i) {
	nodeIDsinTRS.push_back(nodesInFace[i]);
      }
    }

    sort(nodeIDsinTRS.begin(), nodeIDsinTRS.end());
    vector<CFuint> uniqueNodeList;
    uniqueNodeList.reserve(nodeIDsinTRS.size());
    unique_copy(nodeIDsinTRS.begin(),
		nodeIDsinTRS.end(),
		back_inserter(uniqueNodeList));

    const CFuint nbNodes = uniqueNodeList.size();
    cout << "unique TRS nbNodes = " << nbNodes << " , all TRS nbNodes = " << nodeIDsinTRS.size() << endl;
    SafePtr<RealMatrix> nodes = _nodes.getData();

    // map the global node ID to the local ID (in the current TRS)
    CFMap<CFuint,CFuint> mapNodeIDToLocalID(nbNodes);

    fout << "ZONE N=" << nbNodes << ", E= " << nbFaces <<
      ", F=FEPOINT, ET=QUADRILATERAL" << endl;

    for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
      const CFuint nodeID = uniqueNodeList[iNode];
      mapNodeIDToLocalID.insert(nodeID, iNode);
      fout << setw(16) << fixed << setprecision(12)
	   << (*nodes)(nodeID,XX) << " "
	   << (*nodes)(nodeID,YY) << " "
	   << (*nodes)(nodeID,ZZ) << endl;
    }
    mapNodeIDToLocalID.sortKeys();

    for (MapIt f = faces.first; f != faces.second; ++f) {
      int* nodesInFace = f->second.getNodesID();
      for (CFuint i =0; i < 4; ++i) {
        fout << mapNodeIDToLocalID.find(nodesInFace[i]) + 1 << " ";
      }
      fout << endl;
    }

    // reset the list of nodes in TRS
    nodeIDsinTRS.clear();
  }
  cout << endl;

  cout << "totalNbFaces written = " << totalNbFaces << endl;

  fout.close();

  //  cout <<"######### start uniqueness check for boundary faces ########" << endl;
  //   // uniqueness check for boundary faces (highly inefficient implementation!!)
  //   const CFuint nF = _bface.getNbFaces();
  //   for (CFuint iFace = 0; iFace < nF; ++iFace) {
  //     FaceData::Itr faceItr = _bface.getStartData(iFace);
  //     int* fnodes = faceItr.getNodesID();
  //     CFuint countMatchingFaces = 0;
  //     for (CFuint jFace = 0; jFace < nF; ++jFace) {
  //       if (jFace != iFace) {
  // 	FaceData::Itr f = _bface.getStartData(jFace);
  // 	int* otherNodes = f.getNodesID();
  // 	CFuint countMatch = 0;
  // 	// order of nodes doesn't matter here
  // 	for (CFuint i = 0; i < 4; ++i) {
  // 	  for (CFuint j = 0; j < 4; ++j) {
  // 	    if (otherNodes[i] == fnodes[j]) {
  // 	      countMatch++;
  // 	      break;
  // 	    }
  // 	  }
  // 	}

  // 	if (countMatch == 4) {
  // 	  countMatchingFaces++;
  // 	}
  //       }
  //     }

  //     if(countMatchingFaces > 0) {
  //       cout << "found NOT UNIQUE face in BC " << faceItr.getBcID()
  // 	   << ", neighbor elemID = " << faceItr.getNeighborElemID()<< endl;
  //     }
  //   }
  //   cout <<"######### end uniqueness check for boundary faces ########" << endl;

  cout << "ConvertGridProMesh::writeBoundaryTecplotFile() => end" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertGridProMesh::readTopologyFile()
{
  // this routine reads the topology file where all the comments have been
  // already stripped out
  std::string topoFile = _fileName + ".topo";
  cout << "reading file = " <<  topoFile<< endl;
  ifstream fin(topoFile.c_str());
  ofstream fout("out.topo");

  std::string tmp;
  fin >> tmp; // \startopo
  fout << tmp << endl;
  CFVec<CFuint> header(3);
  fin >> header >> _nbBlocks;
  fout << header << " " << _nbBlocks << endl;
  _blocks.resize(_nbBlocks);

  const CFuint nbBlockFaces = 6;
  CFuint blockID = 0;
  CFuint nbSubFaces = 0;
  CFuint subFaceID = 0;
// unused //  CFuint bcType = 0;
  CFuint faceID = 0;
  CFVec<CFuint> idxHost(9);
  CFVec<CFuint> idxDonor(7);

  for (CFuint ib = 0; ib < _nbBlocks; ++ib) {
    fin >> blockID;
    fout << "blockID = " << blockID << endl;
    if (blockID != ib+1) {
      cout << "ERROR: blockID != ib+1"<< endl;
      cout << "blockID = " << blockID << endl;
      abort();
    }
    fin >> tmp; // nx*ny*nz
    fout << "block size = " << tmp << endl;

    // block face fields
    for (CFuint iFace = 0; iFace < nbBlockFaces; ++iFace) {
      fin >> faceID >> tmp >> nbSubFaces;
      if (faceID != iFace+1) {
	cout << "ERROR: faceID != iFace+1"<< endl;
	cout << "block faceID = " << faceID << endl;
	abort();
      }
      if (nbSubFaces != 1) {
	cout << "ERROR: nbSubFaces != 1"<< endl;
	cout << "nbSubFaces = " << nbSubFaces << endl;
	abort();
      }
      fout << "block faceID = " << faceID
	   << ", size face = " << tmp
	   << ", nb sub-faces = " << nbSubFaces << endl;

      fin >> subFaceID >> _blocks[ib].face[iFace].bcID;
      if (subFaceID != 1) {
	cout << "ERROR: subFaceID != 1"<< endl;
	cout << "subFaceID = " << subFaceID << endl;
	abort();
      }
      fout << subFaceID << " "
	   <<  _blocks[ib].face[iFace].bcID << endl;

      fin >> idxHost;
      fout << "host indices = " << idxHost << endl;

      if (_blocks[ib].face[iFace].bcID >= 0) {
	fin >> _blocks[ib].face[iFace].neighBlock
	    >> _blocks[ib].face[iFace].neighFace
	    >> idxDonor;

	fout << "donor indices = "
	     << _blocks[ib].face[iFace].neighBlock << " "
	     << _blocks[ib].face[iFace].neighFace << " "
	     << idxDonor << endl;
      }
    }
  }

  fin >> tmp;
  fout << tmp << " ";
  fin >> tmp;
  fout << tmp << endl;
  fin.close();
}

//////////////////////////////////////////////////////////////////////////////

} // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
