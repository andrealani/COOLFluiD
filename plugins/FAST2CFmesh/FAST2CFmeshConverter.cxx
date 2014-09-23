// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>
#include <set>

#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Common/Stopwatch.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/BadValueException.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"
#include "FAST2CFmesh/FAST2CFmeshConverter.hh"
#include "FAST2CFmesh/FAST2CFmesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FAST2CFmesh {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FAST2CFmeshConverter,
          MeshFormatConverter,
          FAST2CFmeshModule,
          1>
fast2CFmeshConverterProvider("FAST2CFmesh");

//////////////////////////////////////////////////////////////////////////////

FAST2CFmeshConverter::FAST2CFmeshConverter (const std::string& name)
: MeshFormatConverter(name),
  m_dimension(3),
  m_nbVars(0),
  m_nbCells(0),
  m_nbUpdatableNodes(0),
  m_nbFaces(0),
  m_nbPatches(0),
  m_coordinate(CFNULL),
  m_fcon(CFNULL),
  m_bcID(),
  m_mapBcToFaceIDs(),
  m_econ(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

FAST2CFmeshConverter::~FAST2CFmeshConverter()
{
  deletePtr(m_coordinate);
  deletePtr(m_fcon);
  deletePtr(m_econ);
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::readFAST(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
  path meshFile = change_extension(filepath, getOriginExtension());

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle =
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(meshFile);

  // read general information
  fin >> m_nbUpdatableNodes >> m_nbFaces >> m_nbCells;

  // read node coordinates
  m_coordinate = new Table<CFreal>(m_nbUpdatableNodes, m_dimension);
  for (CFuint i = 0; i < m_dimension; ++i) {
    for (CFuint j = 0; j < m_nbUpdatableNodes; ++j) {
      fin >> (*m_coordinate)(j,i);
    }
  }

  // boundary faces connectivity
  m_fcon = new Table<CFuint>(m_nbFaces, m_dimension);
  for (CFuint i = 0; i < m_nbFaces; ++i) {
    for (CFuint j = 0; j < m_dimension; ++j) {
      fin >> (*m_fcon)(i,j);
    }
  }

  // boundary flag
  set<int> uniqueBcIDs;
  vector<int> bcID(m_nbFaces);
  for (CFuint i = 0; i < m_nbFaces; ++i) {
    fin >> bcID[i];
    uniqueBcIDs.insert(bcID[i]);
  }

  m_nbPatches = uniqueBcIDs.size();
  m_bcID.resize(m_nbPatches);

  CFout << "Nb of patches detected = " << m_nbPatches << "\n";

  CFuint iBC = 0;
  set<int>::const_iterator it;
  for (it = uniqueBcIDs.begin(); it != uniqueBcIDs.end(); ++it, ++iBC) {
    m_bcID[iBC] = (*it);
  }

  CFout << "boundary IDs = ";
  for (CFuint i = 0; i < m_bcID.size(); ++i) {
    CFout << m_bcID[i] << " ";
  }
  CFout << "\n";

  m_mapBcToFaceIDs.reserve(m_nbPatches*m_nbFaces);
  for (CFuint i = 0; i < m_nbFaces; ++i) {
    m_mapBcToFaceIDs.insert(bcID[i],i);
  }
  m_mapBcToFaceIDs.sortKeys();

  // element connectivity
  const CFuint nbNodesInElem = m_dimension + 1;
  m_econ = new Table<CFuint>(m_nbCells,  nbNodesInElem);
  for (CFuint i = 0; i < m_nbCells; ++i) {
    for (CFuint j = 0; j <  nbNodesInElem; ++j) {
      fin >> (*m_econ)(i,j);
    }
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeFAST(const boost::filesystem::path& filepath)
{
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!NB_TRSs " << m_nbPatches << "\n";

  typedef CFMultiMap<int,CFuint>::MapIterator MapIterator;

  // only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < m_nbPatches; ++iTRS) {
    const std::string nameTRS = "Side" + Common::StringOps::to_str < CFuint > ( iTRS );
    fout << "!TRS_NAME " << nameTRS << "\n";
    fout << "!NB_TRs "<< 1 << "\n";
    fout << "!NB_GEOM_ENTS ";
    
    bool found = false;
    pair<MapIterator, MapIterator>  mit = m_mapBcToFaceIDs.find(m_bcID[iTRS], found);
    cf_assert(found);
    
    CFuint nbTRFaces = 0;
    for (MapIterator it = mit.first; it != mit.second; ++it) {
      nbTRFaces++;
    }
    fout << nbTRFaces << "\n";

    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";

    for (MapIterator it = mit.first; it != mit.second; ++it) {
      const CFuint nbNodesPerFace = m_dimension;
      const CFuint nbStatesPerFace = m_dimension;
      const CFuint faceID = it->second;

      fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
      for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
  fout << (*m_fcon)(faceID, iNode) -1<< " " ;
      }
      for (CFuint iState = 0; iState < nbStatesPerFace; ++iState) {
  fout << (*m_fcon)(faceID, iState) -1<< " " ;
      }
      fout << "\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeDiscontinuousTrsData(ofstream& fout)
{
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::readFiles(const boost::filesystem::path& filepath)
{
  readFAST(filepath);
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeContinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  fout << "!NB_NODES "
       << m_nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  fout << "!NB_STATES "
       << m_nbUpdatableNodes
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << m_nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << 1 << "\n";

  /// @todo only first order for now
  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER1 << "\n";

  const CFuint nbNodesInElem = m_dimension + 1;

  fout << "!ELEM_TYPES ";
  fout << MapGeoEnt::identifyGeoEnt(nbNodesInElem, CFPolyOrder::ORDER1, m_dimension) << "\n";

  fout << "!NB_ELEM_PER_TYPE ";
  fout << m_nbCells << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  fout << nbNodesInElem << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  fout << nbNodesInElem << "\n";

  fout << "!LIST_ELEM " << "\n";
  for (CFuint i = 0; i < m_nbCells; ++i) {
    for (CFuint j = 0; j < nbNodesInElem; ++j) {
      fout << (*m_econ)(i,j) -1<<" " ;
    }
    for (CFuint j = 0; j < nbNodesInElem; ++j) {
      fout << (*m_econ)(i,j) -1<<" " ;
    }
    fout << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << 0 << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeNodes(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_NODE " << "\n";
  for (CFuint k = 0; k < m_nbUpdatableNodes; ++k) {
    for (CFuint j = 0; j < m_dimension; ++j)
    {
      fout.precision(14);
      fout.setf(ios::scientific,ios::floatfield);
      fout << (*m_coordinate)(k,j) << " ";
    }
    fout << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeDiscontinuousElements(ofstream& fout)
{
}

//////////////////////////////////////////////////////////////////////////////

void FAST2CFmeshConverter::writeDiscontinuousStates(ofstream& fout)
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FAST2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
