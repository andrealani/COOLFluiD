// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConvertStructMesh/ConvertStructMesh.hh"
#include "Common/NoSuchValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////

ConvertStructMesh::ConvertStructMesh(const std::string& fileName) :
  Common::OwnedObject(),
  _fileName(fileName),
  _nbBlocks(0),
  _nbElements(0),
  _imax(0),
  _jmax(0),
  _kmax(0),
  _blocks(),
  _nodes(),
  _elements(),
  _bface(),
  _mapBcIDToFace(),
  _bcListIDs()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvertStructMesh::~ConvertStructMesh()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvertStructMesh::convert(const std::string& meshType)
{
  cout << "ConvertStructMesh::convert() START" << endl;

  readTopologyFile();
  readAndBuildNodes();
  buildElements();
  writeBoundaryTecplotFile();
//  writeTecplotFile();

  if (meshType == "FVM") {
    writeCFmeshFileFVM();
  }
  else if (meshType == "FVMSplit") {
    writeCFmeshFileFVMSplit();
  }
  else if (meshType == "FEM") {
    writeCFmeshFileFEM();
  }
  else if (meshType == "FEMSplit") {
    writeCFmeshFileFEMSplit();
  }
  else {
    cout << "WATCH OUT: ConvertStructMesh::convert() => meshType not specified" << endl;
    writeCFmeshFileFVM();
  }

  cout << "ConvertStructMesh::convert() END" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertStructMesh::writeTecplotFile()
{
  cout << "ConvertStructMesh::writeTecplotFile() => start" << endl;

  std::string fileName = getGlobalMeshFile() + ".plt";
  ofstream fout(fileName.c_str());

  SafePtr<RealMatrix> nodes = _nodes.getData();
  const CFuint nbNodes = nodes->nbRows();

  fout << "TITLE =  Unstructured grid data" << endl;

  if (_kmax > 0) {
    fout << "VARIABLES  =  \"x0\" \"x1\" \"x2\" " << endl;
    fout << "ZONE N=" << nbNodes << ", E= " << _nbElements <<
      ", F=FEPOINT, ET=BRICK" << endl;
  }
  else {
    fout << "VARIABLES  =  \"x0\" \"x1\" " << endl;
    fout << "ZONE N=" << nbNodes << ", E= " << _nbElements <<
      ", F=FEPOINT, ET=QUADRILATERAL" << endl;
  }

  const CFuint dim = (_kmax > 0) ? 3 : 2;
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    for (CFuint in = 0; in < dim; ++in) {
      //   fout << setw(16) << fixed << setprecision(12) << (*nodes)(iNode,in) << " ";
      fout.precision(13);
      fout.setf(ios_base::scientific, ios_base::floatfield);
      fout << (*nodes)(iNode,in) << " ";
    }
    fout <<endl;
  }

  const CFuint nbNodesInElem = (_kmax > 0) ? 8 : 4;
  for (CFuint iElem = 0; iElem < _nbElements; ++iElem) {
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,in) + 1 << " ";
    }
    fout <<endl;
  }

  fout.close();

  cout << "ConvertStructMesh::writeTecplotFile() => end" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertStructMesh::writeCFmeshFileFEM()
{
  cout << "ConvertStructMesh::writeCFmeshFileFEM() => start" << endl;

  std::string fileName = getGlobalMeshFile() + ".FEM.CFmesh";
  ofstream fout(fileName.c_str());

  SafePtr<RealMatrix> nodes = _nodes.getData();
  const CFuint nbNodes = nodes->nbRows();
  const CFuint dim = (_kmax > 0) ? DIM_3D : DIM_2D;
  const CFuint nbNodesInElem = (dim == DIM_3D) ? 8 : 4;

  fout << "!COOLFLUID_VERSION "    << Environment::CFEnv::getInstance().getCFVersion() << "\n";
  // this can fail if there are problems with SVN
  // fout << "!COOLFLUID_SVNVERSION " << Environment::CFEnv::getInstance().getSvnVersion() << "\n";
  fout << "!CFMESH_FORMAT_VERSION 1.2\n";
  fout << "!NB_DIM " << dim << endl;
  fout << "!NB_EQ " << 1 << endl;
  fout << "!NB_NODES " << nbNodes << " " << 0 <<  endl;
  fout << "!NB_STATES " << nbNodes << " " << 0 << endl;
  fout << "!NB_ELEM " << _nbElements << endl;
  fout << "!NB_ELEM_TYPES " << 1 << endl;
  fout << "!GEOM_POLYORDER " << 1 << endl;
  fout << "!SOL_POLYORDER " << 1 << endl;
  if (dim == DIM_3D) {
    fout << "!ELEM_TYPES Hexa" << endl;
  }
  else {
    fout << "!ELEM_TYPES Quad" << endl;
  }

  fout << "!NB_ELEM_PER_TYPE "<< _nbElements << endl;
  fout << "!NB_NODES_PER_TYPE "<< nbNodesInElem << endl;
  fout << "!NB_STATES_PER_TYPE "<< nbNodesInElem << endl;
  fout << "!LIST_ELEM" << endl;
  for (CFuint iElem = 0; iElem < _nbElements; ++iElem) {
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,in) << " ";
      cf_assert(_elements(iElem,in) < nbNodes);
    }
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,in) << " ";
      cf_assert(_elements(iElem,in) < nbNodes);
    }
    fout << endl;
  }

  typedef CFMultiMap<int,FaceData::Itr>::MapIterator MapIt;
  const CFuint nbNodesInFace = (dim == DIM_3D) ? 4 : 2;

  fout << "!NB_TRSs " << _bcListIDs.size() << endl;
  CFuint totalNbFaces = 0;
  for (CFuint iTRS = 0; iTRS < _bcListIDs.size(); ++iTRS) {
    fout << "!TRS_NAME " << "Side" << iTRS << endl;
    fout << "!NB_TRs " << 1 << endl;

    bool fo = false;
    pair<MapIt,MapIt> faces = _mapBcIDToFace.find(_bcListIDs[iTRS], fo);
    cf_assert(fo);
    
    CFuint nbFaces = 0;
    for (MapIt f = faces.first; f != faces.second; ++f) {
      nbFaces++;
    }

    fout << "!NB_GEOM_ENTS " << nbFaces << endl;
    fout << "!GEOM_TYPE Face" << endl;
    fout << "!LIST_GEOM_ENT" << endl;

    const CFuint nbStatesInFace = nbNodesInFace;
    for (MapIt f = faces.first; f != faces.second; ++f, ++totalNbFaces) {
      fout << nbNodesInFace << " " << nbStatesInFace << " ";
      int* nodesInFace = f->second.getNodesID();
      for (CFuint i = 0; i < nbNodesInFace; ++i) {
        fout << nodesInFace[i] << " ";
	cf_assert( (CFuint) nodesInFace[i] < nbNodes);
      }

      for (CFuint i = 0; i < nbNodesInFace; ++i) {
        fout << nodesInFace[i] << " ";
	cf_assert( (CFuint) nodesInFace[i] < nbNodes);
      }
      fout << endl;
      cf_assert(static_cast<CFuint>(f->second.getNeighborElemID()) < _nbElements);
    }
  }

  fout << "!LIST_NODE" << endl;
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    for (CFuint in = 0; in < dim; ++in) {
      //   fout << setw(16) << fixed << setprecision(12) << (*nodes)(iNode,in) << " ";
      fout.precision(13);
      fout.setf(ios_base::scientific, ios_base::floatfield);
      fout << (*nodes)(iNode,in) << " ";
    }
    fout <<endl;
  }
  fout << "!LIST_STATE " << 0 << endl;
  fout << "!END" << endl;

  fout.close();

  cout << "total nb faces = " << totalNbFaces << endl;
  cout << "ConvertStructMesh::writeCFmeshFileFEM() => end" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertStructMesh::writeCFmeshFileFEMSplit()
{
  cout << "ConvertStructMesh::writeCFmeshFileFEMSplit() => start" << endl;

  std::string fileName = getGlobalMeshFile() + ".FEM.Split.CFmesh";
  ofstream fout(fileName.c_str());

  SafePtr<RealMatrix> nodes = _nodes.getData();
  const CFuint nbNodes = nodes->nbRows();
  const CFuint dim = (_kmax > 0) ? DIM_3D : DIM_2D;

  if (dim != DIM_2D) {
    throw NoSuchValueException (FromHere(),"ConvertStructMesh::writeCFmeshFileFEMSplit()=> only 2D available");
  }

  const CFuint nbTriagElements = _nbElements*2;
  const CFuint nbNodesInElem = 3;

  fout << "!COOLFLUID_VERSION "    << Environment::CFEnv::getInstance().getCFVersion() << "\n";
  // this can fail if there are problems with SVN
  // fout << "!COOLFLUID_SVNVERSION " << Environment::CFEnv::getInstance().getSvnVersion() << "\n";
  fout << "!CFMESH_FORMAT_VERSION 1.2\n";
  fout << "!NB_DIM " << dim << endl;

  vector<CFuint> tID1(3);
  tID1[0] = 0; tID1[1] = 1;  tID1[2] = 2;

  vector<CFuint> tID2(3);
  tID2[0] = 0; tID2[1] = 2;  tID2[2] = 3;

  fout << "!NB_EQ " << 1 << endl;
  fout << "!NB_NODES " << nbNodes << " " << 0 <<  endl;
  fout << "!NB_STATES " << nbNodes << " " << 0 << endl;
  fout << "!NB_ELEM " << nbTriagElements << endl;
  fout << "!NB_ELEM_TYPES " << 1 << endl;
  fout << "!GEOM_POLYORDER " << 1 << endl;
  fout << "!SOL_POLYORDER " << 1 << endl;
  fout << "!ELEM_TYPES Triag" << endl;
  fout << "!NB_ELEM_PER_TYPE "<< nbTriagElements << endl;
  fout << "!NB_NODES_PER_TYPE "<< nbNodesInElem << endl;
  fout << "!NB_STATES_PER_TYPE "<< nbNodesInElem << endl;
  fout << "!LIST_ELEM" << endl;
  for (CFuint iElem = 0; iElem < _nbElements; ++iElem) {
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,tID1[in]) << " ";
      cf_assert(_elements(iElem,tID1[in]) < nbNodes);
    }
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,tID1[in]) << " ";
      cf_assert(_elements(iElem,tID1[in]) < nbNodes);
    }
    fout << endl;

    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,tID2[in]) << " ";
      cf_assert(_elements(iElem,tID2[in]) < nbNodes);
    }
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,tID2[in]) << " ";
      cf_assert(_elements(iElem,tID2[in]) < nbNodes);
    }
    fout << endl;
  }

  typedef CFMultiMap<int,FaceData::Itr>::MapIterator MapIt;
  const CFuint nbNodesInFace = 2;

  fout << "!NB_TRSs " << _bcListIDs.size() << endl;
  CFuint totalNbFaces = 0;
  for (CFuint iTRS = 0; iTRS < _bcListIDs.size(); ++iTRS) {
    fout << "!TRS_NAME " << "Side" << iTRS << endl;
    fout << "!NB_TRs " << 1 << endl;
    
    bool fo = false;
    pair<MapIt,MapIt> faces = _mapBcIDToFace.find(_bcListIDs[iTRS], fo);
    cf_assert(fo);
    
    CFuint nbFaces = 0;
    for (MapIt f = faces.first; f != faces.second; ++f) {
      nbFaces++;
    }

    fout << "!NB_GEOM_ENTS " << nbFaces << endl;
    fout << "!GEOM_TYPE Face" << endl;
    fout << "!LIST_GEOM_ENT" << endl;

    const CFuint nbStatesInFace = nbNodesInFace;
    for (MapIt f = faces.first; f != faces.second; ++f, ++totalNbFaces) {
      fout << nbNodesInFace << " " << nbStatesInFace << " ";
      int* nodesInFace = f->second.getNodesID();
      for (CFuint i = 0; i < nbNodesInFace; ++i) {
        fout << nodesInFace[i] << " ";
	cf_assert( (CFuint) nodesInFace[i] < nbNodes);
      }

      for (CFuint i = 0; i < nbNodesInFace; ++i) {
        fout << nodesInFace[i] << " ";
	cf_assert( (CFuint) nodesInFace[i] < nbNodes);
      }
      fout << endl;
    }
  }

  fout << "!LIST_NODE" << endl;
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    for (CFuint in = 0; in < dim; ++in) {
      //   fout << setw(16) << fixed << setprecision(12) << (*nodes)(iNode,in) << " ";
      fout.precision(13);
      fout.setf(ios_base::scientific, ios_base::floatfield);
      fout << (*nodes)(iNode,in) << " ";
    }
    fout <<endl;
  }
  fout << "!LIST_STATE " << 0 << endl;
  fout << "!END" << endl;

  fout.close();

  cout << "total nb faces = " << totalNbFaces << endl;
  cout << "ConvertStructMesh::writeCFmeshFileFEMSplit() => end" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertStructMesh::writeCFmeshFileFVM()
{
  cout << "ConvertStructMesh::writeCFmeshFileFVM() => start" << endl;

  std::string fileName = getGlobalMeshFile() + ".FVM.CFmesh";
  ofstream fout(fileName.c_str());

  SafePtr<RealMatrix> nodes = _nodes.getData();
  const CFuint nbNodes = nodes->nbRows();
  const CFuint dim = (_kmax > 0) ? DIM_3D : DIM_2D;
  const CFuint nbNodesInElem = (dim == DIM_3D) ? 8 : 4;

  fout << "!COOLFLUID_VERSION "    << Environment::CFEnv::getInstance().getCFVersion() << "\n";
  // this can fail if there are problems with SVN
  // fout << "!COOLFLUID_SVNVERSION " << Environment::CFEnv::getInstance().getSvnVersion() << "\n";
  fout << "!CFMESH_FORMAT_VERSION 1.2\n";
  fout << "!NB_DIM " << dim << endl;
  fout << "!NB_EQ " << 1 << endl;
  fout << "!NB_NODES " << nbNodes << " " << 0 <<  endl;
  fout << "!NB_STATES " << _nbElements << " " << 0 << endl;
  fout << "!NB_ELEM " << _nbElements << endl;
  fout << "!NB_ELEM_TYPES " << 1 << endl;
  fout << "!GEOM_POLYORDER " << 1 << endl;
  fout << "!SOL_POLYORDER " << 0 << endl;
  if (dim == DIM_3D) {
    fout << "!ELEM_TYPES Hexa" << endl;
  }
  else {
    fout << "!ELEM_TYPES Quad" << endl;
  }

  fout << "!NB_ELEM_PER_TYPE "<< _nbElements << endl;
  fout << "!NB_NODES_PER_TYPE "<< nbNodesInElem << endl;
  fout << "!NB_STATES_PER_TYPE "<< 1 << endl;
  fout << "!LIST_ELEM" << endl;
  for (CFuint iElem = 0; iElem < _nbElements; ++iElem) {
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,in) << " ";
      cf_assert(_elements(iElem,in) < nbNodes);
    }
    fout << iElem << endl;
  }

  typedef CFMultiMap<int,FaceData::Itr>::MapIterator MapIt;
  const CFuint nbNodesInFace = (dim == DIM_3D) ? 4 : 2;

  fout << "!NB_TRSs " << _bcListIDs.size() << endl;
  CFuint totalNbFaces = 0;
  for (CFuint iTRS = 0; iTRS < _bcListIDs.size(); ++iTRS) {
    fout << "!TRS_NAME " << "Side" << iTRS << endl;
    fout << "!NB_TRs " << 1 << endl;
    
    bool fo = false;
    pair<MapIt,MapIt> faces = _mapBcIDToFace.find(_bcListIDs[iTRS], fo);
    cf_assert(fo);
    
    CFuint nbFaces = 0;
    for (MapIt f = faces.first; f != faces.second; ++f) {
      nbFaces++;
    }

    fout << "!NB_GEOM_ENTS " << nbFaces << endl;
    fout << "!GEOM_TYPE Face" << endl;
    fout << "!LIST_GEOM_ENT" << endl;

    const CFuint nbStatesInFace = 1;
    for (MapIt f = faces.first; f != faces.second; ++f, ++totalNbFaces) {
      fout << nbNodesInFace << " " << nbStatesInFace << " ";
      int* nodesInFace = f->second.getNodesID();
      for (CFuint i = 0; i < nbNodesInFace; ++i) {
        fout << nodesInFace[i] << " ";
	      cf_assert( (CFuint) nodesInFace[i] < nbNodes);
      }
      fout << f->second.getNeighborElemID() << endl;
      cf_assert(static_cast<CFuint>(f->second.getNeighborElemID()) < _nbElements);
    }
  }

  fout << "!LIST_NODE" << endl;
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    for (CFuint in = 0; in < dim; ++in) {
      //   fout << setw(16) << fixed << setprecision(12) << (*nodes)(iNode,in) << " ";
      fout.precision(13);
      fout.setf(ios_base::scientific, ios_base::floatfield);
      fout << (*nodes)(iNode,in) << " ";
    }
    fout <<endl;
  }
  fout << "!LIST_STATE " << 0 << endl;
  fout << "!END" << endl;

  fout.close();

  cout << "total nb faces = " << totalNbFaces << endl;
  cout << "ConvertStructMesh::writeCFmeshFileFVM() => end" << endl;
}

//////////////////////////////////////////////////////////////////////////////

void ConvertStructMesh::writeCFmeshFileFVMSplit()
{
  cout << "ConvertStructMesh::writeCFmeshFileFVMSplit() => start" << endl;

  std::string fileName = getGlobalMeshFile() + ".FVM.Split.CFmesh";
  ofstream fout(fileName.c_str());

  SafePtr<RealMatrix> nodes = _nodes.getData();
  const CFuint nbNodes = nodes->nbRows();
  const CFuint dim = (_kmax > 0) ? DIM_3D : DIM_2D;

  // here we assume to have DIM_2D
  assert(dim == DIM_2D);
  const CFuint nbNodesInElem = 3;
  const CFuint nbTriagElements = _nbElements*2;

  fout << "!COOLFLUID_VERSION "    << Environment::CFEnv::getInstance().getCFVersion() << "\n";
  // this can fail if there are problems with SVN
  // fout << "!COOLFLUID_SVNVERSION " << Environment::CFEnv::getInstance().getSvnVersion() << "\n";
  fout << "!CFMESH_FORMAT_VERSION 1.2\n";
  fout << "!NB_DIM " << dim << endl;

  vector<CFuint> tID1(3);
  tID1[0] = 0; tID1[1] = 1;  tID1[2] = 2;

  vector<CFuint> tID2(3);
  tID2[0] = 0; tID2[1] = 2;  tID2[2] = 3;

  fout << "!NB_EQ " << 1 << endl;
  fout << "!NB_NODES " << nbNodes << " " << 0 <<  endl;
  fout << "!NB_STATES " << nbTriagElements << " " << 0 << endl;
  fout << "!NB_ELEM " << nbTriagElements << endl;
  fout << "!NB_ELEM_TYPES " << 1 << endl;
  fout << "!GEOM_POLYORDER " << 1 << endl;
  fout << "!SOL_POLYORDER " << 0 << endl;
  fout << "!ELEM_TYPES Triag" << endl;
  fout << "!NB_ELEM_PER_TYPE "<< nbTriagElements << endl;
  fout << "!NB_NODES_PER_TYPE "<< nbNodesInElem << endl;
  fout << "!NB_STATES_PER_TYPE "<< 1 << endl;
  fout << "!LIST_ELEM" << endl;

  CFuint countElem = 0;
  for (CFuint iElem = 0; iElem < _nbElements; ++iElem) {
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,tID1[in]) << " ";
      cf_assert(_elements(iElem,tID1[in]) < nbNodes);
    }
    fout << countElem++ << endl;

    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      fout << _elements(iElem,tID2[in]) << " ";
      cf_assert(_elements(iElem,tID2[in]) < nbNodes);
    }
    fout << countElem++ << endl;
  }

  typedef CFMultiMap<int,FaceData::Itr>::MapIterator MapIt;
  const CFuint nbNodesInFace = 2;

  fout << "!NB_TRSs " << _bcListIDs.size() << endl;
  CFuint totalNbFaces = 0;
  for (CFuint iTRS = 0; iTRS < _bcListIDs.size(); ++iTRS) {
    fout << "!TRS_NAME " << "Side" << iTRS << endl;
    fout << "!NB_TRs " << 1 << endl;

    bool fo = false;
    pair<MapIt,MapIt> faces = _mapBcIDToFace.find(_bcListIDs[iTRS],fo);
    cf_assert(fo);
    
    CFuint nbFaces = 0;
    for (MapIt f = faces.first; f != faces.second; ++f) {
      nbFaces++;
    }

    fout << "!NB_GEOM_ENTS " << nbFaces << endl;
    fout << "!GEOM_TYPE Face" << endl;
    fout << "!LIST_GEOM_ENT" << endl;

    const CFuint nbStatesInFace = 1;
    for (MapIt f = faces.first; f != faces.second; ++f, ++totalNbFaces) {
      fout << nbNodesInFace << " " << nbStatesInFace << " ";
      int* nodesInFace = f->second.getNodesID();
      for (CFuint i = 0; i < nbNodesInFace; ++i) {
        fout << nodesInFace[i] << " ";
	cf_assert( (CFuint) nodesInFace[i] < nbNodes);
      }

      CFuint neighElemID = f->second.getNeighborElemID()*2;
      CFuint count = 0;
      for (CFuint in = 0; in < nbNodesInElem; ++in) {
	if ( (int) _elements(f->second.getNeighborElemID(),tID1[in]) == nodesInFace[0]) {
	  count++;
	}
      }
      for (CFuint in = 0; in < nbNodesInElem; ++in) {
	if ( (int) _elements(f->second.getNeighborElemID(),tID1[in]) == nodesInFace[1]) {
	  count++;
	}
      }
      assert(count <= 2);
      if (count == 2) {
	fout << neighElemID << endl;
      }
      else {
	count = 0;
	neighElemID += 1;

	for (CFuint in = 0; in < nbNodesInElem; ++in) {
	  if ( (int) _elements(f->second.getNeighborElemID(),tID2[in]) == nodesInFace[0]) {
	    count++;
	  }
	}
	for (CFuint in = 0; in < nbNodesInElem; ++in) {
	  if ( (int) _elements(f->second.getNeighborElemID(),tID2[in]) == nodesInFace[1]) {
	    count++;
	  }
	}
	assert(count <= 2);
	if (count == 2) {
	  fout << neighElemID << endl;
	}
	else {
	  cout << "neighbor cell NOT FOUND !!" << endl;
	  abort();
	}
      }
      cf_assert(neighElemID < nbTriagElements);
    }
  }

  fout << "!LIST_NODE" << endl;
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    for (CFuint in = 0; in < dim; ++in) {
      //   fout << setw(16) << fixed << setprecision(12) << (*nodes)(iNode,in) << " ";
      fout.precision(13);
      fout.setf(ios_base::scientific, ios_base::floatfield);
      fout << (*nodes)(iNode,in) << " ";
    }
    fout <<endl;
  }
  fout << "!LIST_STATE " << 0 << endl;
  fout << "!END" << endl;

  fout.close();

  cout << "total nb faces = " << totalNbFaces << endl;
  cout << "ConvertStructMesh::writeCFmeshFileFVMSplit() => end" << endl;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
