// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <iostream>
#include <iomanip>

#include "ConvertStructMesh/ConvertQuadMesh.hh"
#include "MathTools/CFMat.hh"
#include "MathTools/CFVec.hh"
#include "Common/CFMap.hh"
#include "ConvertStructMesh/ConvertStructMeshApp.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ConvertQuadMesh, ConvertStructMesh,
			    ConvertStructMeshAppModule, 1>
convertQuadMeshProv("Quad");

//////////////////////////////////////////////////////////////////////////////

ConvertQuadMesh::ConvertQuadMesh(const std::string& fileName) :
  ConvertStructMesh(fileName)
{
  _elemIDs.resize(2);
  FaceData::setSize(DIM_2D);
  _nodes.setNbBlocks(1);
}

//////////////////////////////////////////////////////////////////////////////

ConvertQuadMesh::~ConvertQuadMesh()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvertQuadMesh::readAndBuildNodes()
{
  cout << "ConvertQuadMesh::readAndBuildNodes() => start"<< endl;

  ifstream fin(getGlobalMeshFile().c_str());

  // nodes per block
  _nodes.setNbBlocks(1);
  _nodes.setDim(DIM_2D);

  vector<CFuint>& nxvec = _nodes.getNxVec();
  vector<CFuint>& nyvec = _nodes.getNyVec();
  nxvec.resize(1);
  nyvec.resize(1);

  vector<vector<bool> > isBnode(1);
 // unused // CFuint nbBNodes = 0;

  // maximum number of nodes in x,y,z directions
  fin >> nxvec[0] >> nyvec[0];

  const CFuint totNbNodes = nxvec[0]*nyvec[0];
  isBnode[0].resize(totNbNodes, false);

  const CFuint nx = nxvec[0];
  const CFuint ny = nyvec[0];
  const CFuint imax = nx-1;
  const CFuint jmax = ny-1;
  _nbElements = imax*jmax;

  SafePtr<RealMatrix> nodes = _nodes.getData();
  nodes->resize(totNbNodes, DIM_2D);

  CFuint iNode = 0;
  for (CFuint j = 0; j < ny; ++j) {
    for (CFuint i = 0; i < nx; ++i, ++iNode) {
      cf_assert(iNode < totNbNodes);
      fin >> (*nodes)(iNode,XX) >> (*nodes)(iNode,YY);

      // gory fix to be made optional
      // if we are on the symmetry axis impose strongly to have y=0
      if (i == 0) {
	// cout << "Y before = " << (*nodes)(iNode,YY) << endl;
	(*nodes)(iNode,YY) = 0.0;
	// cout << "Y after  = " << (*nodes)(iNode,YY) << endl;
      }

    }
  }
  _nodes.assignDefaultIDs(0,totNbNodes);

  cout << "ConvertQuadMesh::readAndBuildNodes() => end"<< endl;
}

/////////////////////////////////////////////////////////////

void ConvertQuadMesh::buildElements()
{
  cout << "ConvertQuadMesh::buildElements() => start"<< endl;

  vector<CFuint>& nxvec = _nodes.getNxVec();
  vector<CFuint>& nyvec = _nodes.getNyVec();
// unused //  vector<CFuint>& nzvec = _nodes.getNzVec();

  // build the elements
  _imax = nxvec[0] - 1;
  _jmax = nyvec[0] - 1;

  const CFuint faceCounter = 2*(_imax + _jmax);
  _elements.resize(_nbElements,4);
  _bface.reserve(faceCounter);
  _mapBcIDToFace.reserve(faceCounter);

  CFuint iElem = 0;
  for (CFuint j = 0; j < _jmax; ++j) {
    for (CFuint i = 0; i < _imax; ++i, ++iElem) {
      cf_assert(iElem < _nbElements);

      _elements(iElem,0) = _nodes.getNodeID(0,i,j);
      _elements(iElem,1) = _nodes.getNodeID(0,i+1,j);
      _elements(iElem,2) = _nodes.getNodeID(0,i+1,j+1);
      _elements(iElem,3) = _nodes.getNodeID(0,i,j+1);

      // if the element is a boundary element
      if ((i == 0 || i == _imax-1) ||
	  (j == 0 || j == _jmax-1)) {
	addBoundaryFaces(iElem, i, j);
      }
    }
  }
  _mapBcIDToFace.sortKeys();

  _bcListIDs.reserve(4);
  for (CFuint i = 0; i < 4; ++i) {
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
  cout << "ConvertQuadMesh::buildElements() => end"<< endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertQuadMesh::addBoundaryFaces(CFuint iElem, CFuint i, CFuint j)
{
  // TRS and quadrilateral topological numbering
  //
  //            2
  //    -----------------
  //    |               |
  //j 3 |               | 1
  //    |               |
  //    -----------------
  //            0
  //            i
  if (j == 0) {
    _elemIDs[0] = _elements(iElem,0);
    _elemIDs[1] = _elements(iElem,1);
    _mapBcIDToFace.insert(0,_bface.addData(0,iElem, 0,_elemIDs));
  }

  if (i == _imax-1) {
    _elemIDs[0] = _elements(iElem,1);
    _elemIDs[1] = _elements(iElem,2);
    _mapBcIDToFace.insert(1,_bface.addData(1,iElem, 1,_elemIDs));
  }

  if (j == _jmax-1) {
    _elemIDs[0] = _elements(iElem,2);
    _elemIDs[1] = _elements(iElem,3);
    _mapBcIDToFace.insert(2,_bface.addData(2,iElem, 2,_elemIDs));
  }

  if (i == 0) {
    _elemIDs[0] = _elements(iElem,3);
    _elemIDs[1] = _elements(iElem,0);
    _mapBcIDToFace.insert(3,_bface.addData(3,iElem, 3,_elemIDs));
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
