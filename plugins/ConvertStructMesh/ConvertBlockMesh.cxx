// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <iostream>
#include <iomanip>

#include <boost/filesystem/path.hpp>

#include "MathTools/CFMat.hh"
#include "MathTools/CFVec.hh"
#include "Common/CFMap.hh"
#include "Environment/ObjectProvider.hh"

#include "ConvertStructMesh/ConvertBlockMesh.hh"
#include "ConvertStructMesh/ConvertStructMeshApp.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ConvertBlockMesh, ConvertStructMesh,
			    ConvertStructMeshAppModule, 1>
convertBlockMeshProv("Block");

//////////////////////////////////////////////////////////////////////////////

ConvertBlockMesh::ConvertBlockMesh(const std::string& fileName) :
  ConvertStructMesh(fileName)
{
  _elemIDs.resize(4);
  FaceData::setSize(DIM_3D);
  _nodes.setNbBlocks(1);
}

//////////////////////////////////////////////////////////////////////////////

ConvertBlockMesh::~ConvertBlockMesh()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvertBlockMesh::readAndBuildNodes()
{
  cout << "ConvertBlockMesh::readAndBuildNodes() => start"<< endl;

  ifstream fin(getGlobalMeshFile().c_str());

  // nodes per block
  _nodes.setNbBlocks(1);
  _nodes.setDim(DIM_3D);

  vector<CFuint>& nxvec = _nodes.getNxVec();
  vector<CFuint>& nyvec = _nodes.getNyVec();
  vector<CFuint>& nzvec = _nodes.getNzVec();

  nxvec.resize(1);
  nyvec.resize(1);
  nzvec.resize(1);

  vector<vector<bool> > isBnode(1);

  // maximum number of nodes in x,y,z directions
  fin >> nxvec[0] >> nyvec[0] >> nzvec[0];

  const CFuint totNbNodes = nxvec[0]*nyvec[0]*nzvec[0];
  isBnode[0].resize(totNbNodes, false);

  const CFuint nx = nxvec[0];
  const CFuint ny = nyvec[0];
  const CFuint nz = nzvec[0];

  const CFuint imax = nx-1;
  const CFuint jmax = ny-1;
  const CFuint kmax = nz-1;

  _nbElements = imax*jmax*kmax;

  SafePtr<RealMatrix> nodes = _nodes.getData();
  nodes->resize(totNbNodes, DIM_3D);

  CFuint iNode = 0;
  for (CFuint k = 0; k < nz; ++k) {
    for (CFuint j = 0; j < ny; ++j) {
      for (CFuint i = 0; i < nx; ++i, ++iNode) {
	cf_assert(iNode < totNbNodes);
	fin >> (*nodes)(iNode,XX) >> (*nodes)(iNode,YY) >> (*nodes)(iNode,ZZ);
      }
    }
  }
  _nodes.assignDefaultIDs(0,totNbNodes);

  cout << "ConvertBlockMesh::readAndBuildNodes() => end"<< endl;
}

/////////////////////////////////////////////////////////////

void ConvertBlockMesh::buildElements()
{
  cout << "ConvertBlockMesh::buildElements() => start"<< endl;

  vector<CFuint>& nxvec = _nodes.getNxVec();
  vector<CFuint>& nyvec = _nodes.getNyVec();
  vector<CFuint>& nzvec = _nodes.getNzVec();

  // build the elements
  _imax = nxvec[0] - 1;
  _jmax = nyvec[0] - 1;
  _kmax = nzvec[0] - 1;

  // total number of boundary faces
  const CFuint faceCounter = 2*((_imax*_jmax) + (_imax*_kmax) + (_jmax*_kmax));
  _elements.resize(_nbElements,8);
  _bface.reserve(faceCounter);
  _mapBcIDToFace.reserve(faceCounter);

  CFuint iElem = 0;
  for (CFuint k = 0; k < _kmax; ++k) {
    for (CFuint j = 0; j < _jmax; ++j) {
      for (CFuint i = 0; i < _imax; ++i, ++iElem) {
	cf_assert(iElem < _nbElements);

	// watch out for the numbering
	_elements(iElem,0) = _nodes.getNodeID(0,i,j,k);
	_elements(iElem,1) = _nodes.getNodeID(0,i+1,j,k);
	_elements(iElem,5) = _nodes.getNodeID(0,i+1,j+1,k);
	_elements(iElem,4) = _nodes.getNodeID(0,i,j+1,k);
	_elements(iElem,3) = _nodes.getNodeID(0,i,j,k+1);
	_elements(iElem,2) = _nodes.getNodeID(0,i+1,j,k+1);
	_elements(iElem,6) = _nodes.getNodeID(0,i+1,j+1,k+1);
	_elements(iElem,7) = _nodes.getNodeID(0,i,j+1,k+1);

	// if the element is a boundary element
	if ((i == 0 || i == _imax-1) ||
	    (j == 0 || j == _jmax-1) ||
	    (k == 0 || k == _kmax-1)) {
	  addBoundaryFaces(iElem, i, j, k);
	}
      }
    }
  }

  _mapBcIDToFace.sortKeys();

  const CFuint nbTopologicalBCFaces = 6;
  _bcListIDs.reserve(nbTopologicalBCFaces);
  for (CFuint i = 0; i < nbTopologicalBCFaces; ++i) {
    _bcListIDs.push_back(i);
  }

  if (_bface.getNbFaces() !=  faceCounter) {
    cout << "ERROR : _bface.getNbFaces() != faceCounter" << endl;
    cout << "_bface.getNbFaces() = " << _bface.getNbFaces() << endl;
    cout << "faceCounter = " << faceCounter << endl;
    abort();
  }
  else {
    cout << "nb bfaces detected = " << _bface.getNbFaces() << endl;
  }

  cout << _nbElements << " elements built" << endl;
  cout << "ConvertBlockMesh::buildElements() => end"<< endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertBlockMesh::addBoundaryFaces(CFuint iElem, CFuint i, CFuint j, CFuint k)
{
  // COOLFLUID numbering for Hexa is used

  if (j == 0) {
    _elemIDs[0] = _elements(iElem,0);
    _elemIDs[1] = _elements(iElem,3);
    _elemIDs[2] = _elements(iElem,2);
    _elemIDs[3] = _elements(iElem,1);
    _mapBcIDToFace.insert(0,_bface.addData(0,iElem, 0,_elemIDs));
  }

  if (j == _jmax-1) {
    _elemIDs[0] = _elements(iElem,4);
    _elemIDs[1] = _elements(iElem,5);
    _elemIDs[2] = _elements(iElem,6);
    _elemIDs[3] = _elements(iElem,7);
    _mapBcIDToFace.insert(1,_bface.addData(0,iElem, 1,_elemIDs));
  }

  if (i == 0) {
    _elemIDs[0] = _elements(iElem,3);
    _elemIDs[1] = _elements(iElem,0);
    _elemIDs[2] = _elements(iElem,4);
    _elemIDs[3] = _elements(iElem,7);
    _mapBcIDToFace.insert(5,_bface.addData(0,iElem, 5,_elemIDs));
  }

  if (i == _imax -1) {
    _elemIDs[0] = _elements(iElem,1);
    _elemIDs[1] = _elements(iElem,2);
    _elemIDs[2] = _elements(iElem,6);
    _elemIDs[3] = _elements(iElem,5);
    _mapBcIDToFace.insert(3,_bface.addData(0,iElem, 3,_elemIDs));
  }

  if (k == 0) {
    _elemIDs[0] = _elements(iElem,0);
    _elemIDs[1] = _elements(iElem,1);
    _elemIDs[2] = _elements(iElem,5);
    _elemIDs[3] = _elements(iElem,4);
    _mapBcIDToFace.insert(2,_bface.addData(0,iElem, 2,_elemIDs));
  }

  if (k == _kmax -1) {
    _elemIDs[0] = _elements(iElem,2);
    _elemIDs[1] = _elements(iElem,3);
    _elemIDs[2] = _elements(iElem,7);
    _elemIDs[3] = _elements(iElem,6);
    _mapBcIDToFace.insert(4,_bface.addData(0,iElem, 4,_elemIDs));
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
