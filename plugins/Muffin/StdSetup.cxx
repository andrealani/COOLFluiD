
#include <boost/progress.hpp>
#include "Common/CFMultiMap.hh"
#include "Common/Table.hh"
#include "Common/NullPointerException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/BC.hh"
#include "Muffin/Loop.hh"
#include "Muffin/StdSetup.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,MuffinData,MuffinModule > cStdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  MuffinCom(name),
  s_nodes("nodes"),                        // socket sinks
  s_states("states"),                      // ...
  s_bstatesneighbors("bStatesNeighbors"),  // ...
  s_nstatesproxy("nstatesProxy"),          // socket sources
  s_faceneighcell("faceNeighCell"),        // ...
  s_mn_volume("NodalVolume"),              // ...
  s_mn_bnormal("NodalNormals"),            // ...
  s_mn_bnarea("NodalArea"),                // ...
  s_mn_walldistance("NodalWallDistance"),  // ...
  s_mn_wallnode("NodalWallNode"),          // ...
  s_mn_priority("NodalPriority")           // ...
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

//std::vector< SafePtr< TopologicalRegionSet > > vtrs_inner = MeshDataStack::getActive()->getFilteredTrsList("inner");
  std::vector< SafePtr< TopologicalRegionSet > > vtrs_boundary = MeshDataStack::getActive()->getFilteredTrsList("boundary");


  d.log("set master loop...");
  {
    if (!d.m_vcomm_loops.size())
      d.err("didn't find a Loop to start from!");
    // set first loop as master
    SafePtr< Loop > l = *d.m_vcomm_loops.begin();
    l->setMaster();
    d.log("master loop: \"" + l->getName() + "\"");
  }
  d.ver("set master loop.");


  d.log("set nodal states proxy...");
  setNStatesProxy();
  d.ver("set nodal states proxy.");


  d.log("set faces neighboring cells...");
  setFaceNeighCell(vtrs_boundary);
  d.ver("set faces neighboring cells.");


  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< std::valarray< State* > > h_bstatesneighbors = s_bstatesneighbors.getDataHandle();
  DataHandle< std::pair< CFuint,CFuint > > h_faceneighcell = s_faceneighcell.getDataHandle();

  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();


  d.log("mark walls...");
  {
    // get TRS names
    std::vector< std::string > walls_names;
    for (CFuint i=0; i<d.m_vcomm_bc.size(); ++i) {
      if (d.m_vcomm_bc[i]->hasTag("Muffin::BCWall")) {
        const std::vector< std::string > wnames = d.m_vcomm_bc[i]->getTrsNames();
        walls_names.insert(walls_names.end(),wnames.begin(),wnames.end());
      }
    }

    // remove duplicate names
    sort(walls_names.begin(),walls_names.end());
    walls_names.erase(
      unique(walls_names.begin(),walls_names.end()),
      walls_names.end() );

    // get TopologicalRegionSet's from names
    d.m_walls.clear();
    for (CFuint i=0; i<walls_names.size(); ++i)
      d.m_walls.push_back(MeshDataStack::getActive()->getTrs(walls_names[i]));
  }
  d.ver("mark walls.");


  d.log("set common variables...");
  nulam     = d.m_kviscosity;
  diffusion = (d.m_diffusion? 1:0);
  buoyancy  = (d.m_buoyancy?  1:0);

  turmod    = ITNULL;
  turmod_ke = false;
  if (     d.m_turbulencemodel=="KE2L") { turmod = ITKE2L; turmod_ke = true; }
  else if (d.m_turbulencemodel=="KELB") { turmod = ITKELB; turmod_ke = true; }
  else if (d.m_turbulencemodel=="KENA") { turmod = ITKENA; turmod_ke = true; }
  else if (d.m_turbulencemodel=="KWHR") { turmod = ITKWHR; }
  else if (d.m_turbulencemodel=="KWLR") { turmod = ITKWLR; }
  else if (d.m_turbulencemodel=="KWPD") { turmod = ITKWPD; }
  else if (d.m_turbulencemodel=="SST")  { turmod = ITKWSS; }
  else if (d.m_turbulencemodel=="BSL")  { turmod = ITKWBS; }
  turmod_walldistance = (
    turmod==ITKE2L || turmod==ITKELB || turmod==ITKENA ||
    turmod==ITKWSS || turmod==ITKWBS ?  1:0 );

  Nvtcell = Ndim+1;
  Nvtfce  = Ndim;
  geo2nodes = MeshDataStack::getActive()->getTrs("InnerCells")->getGeo2NodesConn();
  cf_assert(geo2nodes.isNotNull());
  d.ver("set common variables.");


  d.log("get boundary nodes' neighbors...");
  {
    neighbors.resize(h_nodes.size());
    for (CFuint i=0; i<vtrs_boundary.size(); ++i) {

      SafePtr< std::vector< CFuint > > trsStates = vtrs_boundary[i]->getStatesInTrs();
      std::vector< CFuint >::iterator itd;
      for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
        const CFuint nLocalID = *itd;
        const CFuint nbNeigh = h_bstatesneighbors[nLocalID].size();
        if (!nbNeigh)
          continue;

        // add the neighbors (n) to this node (nLocalID)...
        neighbors[nLocalID].reserve(nbNeigh);
        for (CFuint n=0; n<nbNeigh; ++n)
          neighbors[nLocalID].push_back(
            (int) h_bstatesneighbors[nLocalID][n]->getLocalID() );

        // and remove duplicate neighbors
        sort(neighbors[nLocalID].begin(),neighbors[nLocalID].end());
        neighbors[nLocalID].erase(
          unique( neighbors[nLocalID].begin(), neighbors[nLocalID].end() ),
          neighbors[nLocalID].end() );
      }
    }
  }
  d.ver("get boundary nodes' neighbors.");


  d.m_wall_distance = (turmod_walldistance || d.m_wall_distance);
  if (d.m_wall_distance) {
    d.log("set wall distance...");
    setWallDistance(d.m_walls);
    d.ver("set wall distance.");
  }


  d.log("set boundaries nodal priorities...");
  h_mn_priority.resize(h_nodes.size());
  h_mn_priority = 0;
  for (CFuint i=0; i<d.m_vcomm_bc.size(); ++i) {
    const std::vector< std::string > trs_names = d.m_vcomm_bc[i]->getTrsNames();
    for (CFuint j=0; j<trs_names.size(); ++j) {
      const std::vector< CFuint >& nodes = *MeshDataStack::getActive()->getTrs(trs_names[j])->getNodesInTrs();
      for (CFuint n=0; n<nodes.size(); ++n)
        h_mn_priority[ nodes[n] ] = std::max(
          (CFuint) d.m_vcomm_bc[i]->m_bnpriority,
          h_mn_priority[nodes[n]] );
    }
  }
  d.ver("set boundaries nodal priorities.");


  d.log("set boundaries nodal normals...");
  setBNodeNormals(vtrs_boundary);
  d.ver("set boundary nodal normals.");


  d.log("calculate volume...");
  setNodeVolume();
  d.ver("calculate volume.");


  if (!d.m_restart) {
    d.log("set initial solution...");
    setInitialSolution();
    d.ver("set initial solution.");
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setNStatesProxy()
{
  CFAUTOTRACE;

  // DataHandle to set
  DataHandle< ProxyDofIterator< RealVector >* > h_nstatesproxy = s_nstatesproxy.getDataHandle();

  // get states DataHandle (and its size)
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  const CFuint nbStates = h_states.size();

  //set node to state mapping
  m_nodeIDToStateID.resize(nbStates);
  for (CFuint stateID=0; stateID<nbStates; ++stateID) {
    const CFuint nodeID = h_states[stateID]->getCoordinates().getLocalID();
    cf_assert(nodeID<nbStates);
    m_nodeIDToStateID[nodeID] = stateID;
  }

  h_nstatesproxy.resize(1);
  h_nstatesproxy[0] = new DofDataHandleIterator< RealVector,State,GLOBAL >(h_states, &m_nodeIDToStateID);
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setFaceNeighCell(const std::vector< SafePtr< TopologicalRegionSet > >& btrs)
{
  CFAUTOTRACE;

  // DataHandle to set
  DataHandle<std::pair< CFuint,CFuint > > h_faceneighcell = s_faceneighcell.getDataHandle();

  // get other DataHandle's of interest
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  // get MeshData information
  const TopologicalRegionSet& cells                 = *MeshDataStack::getActive()->getTrs("InnerCells");
  const MapGeoToTrsAndIdx& mapGeoToTrs              = *MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  const std::vector< ElementTypeData >& elementType = *MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType.size();
  const CFuint nbNodes = h_nodes.size();


  std::vector< bool > isBNode(nbNodes,false);  // if node is on the boundary
  CFuint nbBFaces = 0;                         // count boundary faces
  CFMultiMap< CFuint,CFuint > mapBNode2BFace;  // map boundary nodes to faces
  for (CFuint i=0; i<btrs.size(); ++i) {

    const CFuint nbTrsFaces = btrs[i]->getLocalNbGeoEnts();
    nbBFaces += nbTrsFaces;                    // count boundary faces
    for (CFuint iFace = 0; iFace<nbTrsFaces; ++iFace) {
      const CFuint faceID = btrs[i]->getLocalGeoID(iFace);
      for (CFuint iNode = 0; iNode<btrs[i]->getNbNodesInGeo(iFace); ++iNode) {
        const CFuint nodeID = btrs[i]->getNodeID(iFace,iNode);
        cf_assert(nodeID<isBNode.size());
        mapBNode2BFace.insert(nodeID,faceID);  // associate nodeID to faceID
        isBNode[nodeID] = true;                // mark boundary node
      }
    }

  }
  mapBNode2BFace.sortKeys();


  // local connectivity face-node for each element type
  std::vector< Table< CFuint >* > faceNodeElement(nbElemTypes);
  for (CFuint iType = 0; iType<nbElemTypes; ++iType) {
    faceNodeElement[iType] = LocalConnectionData::getInstance().getFaceDofLocal(
       elementType[iType].getGeoShape(),
       static_cast< CFPolyOrder::Type >(elementType[iType].getGeoOrder()),
       NODE,
       CFPolyForm::LAGRANGE );
  }

  // set the face shapes per element type
  std::vector< std::vector< CFGeoShape::Type > > faceShapesPerElemType(nbElemTypes);
  LocalConnectionData::getInstance().setFaceShapesPerElemType(faceShapesPerElemType);

  // number of nodes in a face (preallocated to avoid reallocations)
  std::vector< CFuint > nodesInFace(100);

  // loop over the elements and construct faceIDs
  typedef CFMultiMap< CFuint,CFuint >::MapIterator MapIterator;

  // set h_faceneighcell, looping over the types
  h_faceneighcell.resize(nbBFaces);
  CFuint elemID = 0;
  CFuint countF = 0;
  for (CFuint iType=0; iType<nbElemTypes; ++iType) {
    /// @todo for now all geoents have same geometric and solution polyorder
    const CFuint nbElemFaces = faceShapesPerElemType[iType].size();

    // loop over the elements of this type
    const CFuint nbElemPerType = elementType[iType].getNbElems();

    for (CFuint iElem=0; iElem<nbElemPerType; ++iElem, ++elemID) {
      // loop over the faces in the current element
      for (CFuint iFace=0; iFace<nbElemFaces; ++iFace) {
        const CFuint nbNodesPerFace = faceNodeElement[iType]->nbCols(iFace);

        // construct sets of nodes that make the corresponding face
        // in this element
        CFuint countFaceNodesOnBoundary = 0;
        for (CFuint iNode=0; iNode<nbNodesPerFace; ++iNode) {
          const CFuint localNodeID = (*faceNodeElement[iType])(iFace, iNode);
          const CFuint nodeID = cells.getNodeID(elemID,localNodeID);
          nodesInFace[iNode] = nodeID;
          if (isBNode[nodeID]) {
            ++countFaceNodesOnBoundary;
          }
          else
            break;
        }

        // face is on the boundary if ALL its nodes are on the boundary
        if (countFaceNodesOnBoundary==nbNodesPerFace) {
          // consider the first node belonging to the current face
          // check if you find a face ID shared between all the other
          // nodes building a face
          const CFuint nodeID = nodesInFace[0];
          bool found = false;
	  std::pair< MapIterator,MapIterator > faces = mapBNode2BFace.find(nodeID, found);
	  cf_assert(found);
	  
          // loop over all the faceIDs referenced by the first node to see if
          // all the other nodes reference the same face
          bool faceFound = false;
          CFuint currFaceID = 0;
          for (MapIterator faceInMapItr=faces.first;
               faceInMapItr != faces.second && (faceFound==false);
               ++faceInMapItr) {

            currFaceID = faceInMapItr->second;
            SafePtr< TopologicalRegionSet > faceTrs = mapGeoToTrs.getTrs(currFaceID);

            const CFuint faceIdx = mapGeoToTrs.getIdxInTrs(currFaceID);
            const CFuint nbNodesInCurrFace = faceTrs->getNbNodesInGeo(faceIdx);

            CFuint countNodes = 0;
            // first check if the number of nodes of the face is equal
            if (nbNodesInCurrFace==nbNodesPerFace) {
              for (CFuint iNode=0; iNode<nbNodesPerFace; ++iNode) {
                const CFuint newNodeID = nodesInFace[iNode];
                for (CFuint jNode=0; jNode<nbNodesPerFace; ++jNode) {
                  if (faceTrs->getNodeID(faceIdx,jNode)==newNodeID) {
                    ++countNodes;
                    break;
                  }
                }
              }
            }

            if (countNodes==nbNodesPerFace) {
              // set the neighbor cell ID
              h_faceneighcell[currFaceID].first = elemID;
              // set the local face ID in the corresponding cell
              h_faceneighcell[currFaceID].second = iFace;
              faceFound = true;
              ++countF;
            }
          }

          /*
           * For "corner" cells, it could happen a face whose all nodes are on
           * the boundary but belong to different TRs (possibly even within
           * the same TRS).  In this case the face is NOT a boundary face and
           * so at this point faceFound would remain "false"
           */
          if (!faceFound) {
            /// TODO a more rigorous check could be made to verify
            /// that this exceptional situation is occurring
            CFLogDebugMin("corner inner face with all boundary nodes found\n");
            CFLogDebugMin("face nodes are: ");
            for (CFuint in = 0; in<nbNodesPerFace; ++in) {
              CFLogDebugMin(nodesInFace[in] << " ");
            }
            CFLogDebugMin("\n");
          }
        }
      }
    }
  }

  if (countF!=nbBFaces) {
    throw NoSuchValueException(FromHere(),"not all the face neighbor cell have been detected correctly!!");
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setWallDistance(const std::vector< SafePtr< TopologicalRegionSet > >& walls)
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  const DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const CFuint nbNodes = h_nodes.size();
  const CFuint nbDim = (h_nodes.size()? h_nodes[0]->size() : 0);


  // allocate DataHandle's and initialize
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();
  h_mn_walldistance.resize(nbNodes);
  h_mn_walldistance = 0.;

  DataHandle< CFuint > h_mn_wallnode = s_mn_wallnode.getDataHandle();
  h_mn_wallnode.resize(nbNodes);
  h_mn_wallnode = 0;


  d.log("marking wall nodes...");
  std::vector< bool > nodeiswall(nbNodes,false);
  for (CFuint i=0; i<walls.size(); ++i) {
    std::vector< CFuint >& wallnodes = *walls[i]->getNodesInTrs();
    for (CFuint inb=0; inb<wallnodes.size(); ++inb)
      nodeiswall[wallnodes[inb]] = true;
  }
  d.ver("marking wall nodes.");


  // node to wall distance (minimum distance to each wall node)
  //FIXME make a simple wallnodes vector, it's simpler to parallelize
  d.log("calculating wall distance...");
  if (d.isParallel()) {

    // for each global index
    const CFuint GNnode = MeshDataStack::getActive()->getTotalNodeCount();
    boost::progress_display pbar(GNnode);
    for (CFuint n=0; n<GNnode; ++n, ++pbar) {

      // rank with node sets the coordinates of the node
      int localID = -1;
      double x = -1.e99;
      double y = -1.e99;
      double z = -1.e99;
      {
        for (CFuint o=0; o<nbNodes; ++o)
          if (h_nodes[o]->getGlobalID() == n) {
            x = (*h_nodes[o])[0];
            y = (*h_nodes[o])[1];
            z = (nbDim==2? 0.:(*h_nodes[o])[2]);
            localID = (int) o;
            break;
          }
        GlobalReduceOperation< GRO_MAX >(&x,&x);
        GlobalReduceOperation< GRO_MAX >(&y,&y);
        GlobalReduceOperation< GRO_MAX >(&z,&z);
      }

      const std::string nsp = getMethodData().getNamespace();
      
      // for each wall
      PE::GetPE().setBarrier(nsp);
      double dis = 1.e99;
      int    o = -1;
      for (CFuint i=0; i<walls.size(); ++i) {
        std::vector< CFuint >& wallnodes = *walls[i]->getNodesInTrs();

        // calculate distance to each face node
        for (CFuint inb=0; inb<wallnodes.size(); ++inb) {
          if (localID==(int) wallnodes[inb])  continue;
          const Node &P = (*h_nodes[ wallnodes[inb] ]);
          const double _d = d.distance(x,y,z, P[0],P[1],(nbDim==2?0.:P[2]));
          if (_d<dis) {
            dis = _d;
            o   = wallnodes[inb];
          }
        }

      }

      PE::GetPE().setBarrier(nsp);
      GlobalReduceOperation< GRO_MIN >(&dis,&dis);
      if (localID>0) {
        h_mn_walldistance[localID] = dis;
        h_mn_wallnode[localID] = o;
      }
    }

  }
  else {

    // for each global index
    boost::progress_display pbar(nbNodes);
    for (CFuint n=0; n<nbNodes; ++n, ++pbar) {

      const double x = (*h_nodes[n])[0];
      const double y = (*h_nodes[n])[1];
      const double z = (nbDim==3? (*h_nodes[n])[2]:0.);

      // for each wall
      double dis = 1.e99;
      int    o = -1;
      for (CFuint i=0; i<walls.size(); ++i) {
        const std::vector< CFuint >& wallnodes = *walls[i]->getNodesInTrs();

        // calculate distance to each face node
        for (CFuint inb=0; inb<wallnodes.size(); ++inb) {
          if (n==wallnodes[inb])  continue;
          const Node &P = (*h_nodes[wallnodes[inb]]);
          const double _d = d.distance(x,y,z, P[0],P[1],(nbDim==3? P[2]:0.));
          if (_d<dis) {
            dis = _d;
            o = (int) wallnodes[inb];
          }
        }

      }
      h_mn_walldistance[n] = dis;
      h_mn_wallnode[n] = o;
    }

  }
  d.ver("calculating wall distance.");


  d.log("assigning wall nodes distance...");
  std::vector< int > cornernodes;
  for (CFuint i=0; i<walls.size(); ++i) {
    std::vector< CFuint >& wallnodes = *walls[i]->getNodesInTrs();
    for (CFuint inb=0; inb<wallnodes.size(); ++inb) {
      const int n = wallnodes[inb];
      double dis = 1.e99;
      int    o = -1;

      const double x = (*h_nodes[n])[0];
      const double y = (*h_nodes[n])[1];
      const double z = (nbDim==3? (*h_nodes[n])[2]:0.);

      // use only nodes which are not a wall
      for (CFuint i=0; i<neighbors[n].size(); ++i) {
        const int nn = neighbors[n][i];
        if (nodeiswall[nn])
          continue;

        const Node &P = (*h_nodes[ nn ]);
        const double _d = d.distance(x,y,z, P[0],P[1],(nbDim==3? P[2]:0.));
        if (_d<dis) {
          dis = _d;
          o = nn;
        }

      }
      h_mn_walldistance[n] = dis;
      h_mn_wallnode[n] = o;
      if (o<0)
        cornernodes.push_back(n);
    }
  }

  // assigning corner nodes distance (done after the previous loop because
  // it references calculated distances of the neighbor nodes
  for (CFuint i=0; i<cornernodes.size(); ++i) {
    const int n = cornernodes[i];
    double _d = 1.e99;
    int    o = -1;

    for (CFuint i=0; i<neighbors[n].size(); ++i) {
      const int nn = neighbors[n][i];
      if (n!=nn && h_mn_walldistance[nn]<_d) {
        _d = h_mn_walldistance[nn];
        o  = nn;
      }
    }
    h_mn_walldistance[n] = _d;
    h_mn_wallnode[n] = o;
  }
  d.ver("assigning wall nodes distance.");


#if 0
  d.log("output wall distance...");
  const CFuint nbCells = geo2nodes->nbRows();
  std::ofstream f(d.getFilename("debug-wall-distance.plt").c_str(),std::ios::trunc);
  f << "VARIABLES = x y "<< (nbDim>2? 'z':' ') << " d n" << std::endl;
  f << "ZONE DATAPACKING=POINT ZONETYPE=" << (nbDim>2? "FETETRAHEDRON":"FETRIANGLE") << " N=" << nbNodes << " E=" << nbCells << std::endl;
  for (CFuint n=0; n<nbNodes; ++n)
    f << *h_nodes[n] << ' ' << h_mn_walldistance[n] << ' ' << h_mn_wallnode[n] << std::endl;
  for (CFuint ic=0; ic<nbCells; ++ic) {
    for (int i=0; i<Nvtcell; ++i)
      f << " " << (*geo2nodes)(ic,i)+1;
    f << std::endl;
  }
  f.close();
  d.ver("output wall distance.");
#endif
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setBNodeNormals(const std::vector< SafePtr< TopologicalRegionSet > >& btrs)
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< std::pair< CFuint,CFuint > > h_faceneighcell = s_faceneighcell.getDataHandle();
  const CFuint nbNodes = h_nodes.size();
  const CFuint nbDim = (h_nodes.size()? h_nodes[0]->size() : 0);


  // initialize
  DataHandle< RealVector > h_mn_bnormal = s_mn_bnormal.getDataHandle();
  h_mn_bnormal.resize(nbNodes);
  for (CFuint n=0; n<h_mn_bnormal.size(); ++n)
    h_mn_bnormal[n].resize(nbDim,0.);
  DataHandle< CFreal > h_mn_bnarea = s_mn_bnarea.getDataHandle();
  h_mn_bnarea.resize(nbNodes);
  h_mn_bnarea = 0.;


  struct local_node_struct No_loc[4];
  double vol;
  int inc_min;
  RealVector normal(0.,nbDim);


  // for each boundary face...
  for (CFuint i=0; i<btrs.size(); ++i) {
    const ConnectivityTable< CFuint > faces = *btrs[i]->getGeo2NodesConn();
    for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc) {

      // get inner cell, cell face index and opposite node index
      const CFuint c = h_faceneighcell[ btrs[i]->getLocalGeoID(ifc) ].first;
      const CFuint f = h_faceneighcell[ btrs[i]->getLocalGeoID(ifc) ].second;
      const int o = (nbDim==2? (f==0? 2 : (f==1? 0 : (f==2? 1 :            -1 )))  :
                    (nbDim==3? (f==0? 3 : (f==1? 2 : (f==2? 0 : (f==3? 1 : -1 )))) :
                                                                           -1 ));
      cf_assert_desc("impossible face!",o>=0);

      // calculate normal to face
      d.cellgeom(c,No_loc,&vol,&inc_min);

      // contribute to (boundary) nodal area
      const CFreal fsize = sqrt(No_loc[o].norm2);
      for (int vf=0; vf<Nvtfce; ++vf)
        h_mn_bnarea[ faces(ifc,vf) ] += fsize/(double) Nvtfce;

      // contribute to (boundary) nodal normal
      for (int vf=0; vf<Nvtfce; ++vf) {
        for (CFuint d=0; d<nbDim; ++d)
          normal[d] = No_loc[o].norm[d];
        h_mn_bnormal[ faces(ifc,vf) ] += normal;
      }

    }
  }


  // normalize
  for (CFuint n=0; n<h_mn_bnormal.size(); ++n) {
    if (h_mn_bnormal[n].norm1())
      h_mn_bnormal[n].normalize();
  }


#if 0
  d.log("output boundary normals...");
  std::ofstream f(d.getFilename("debug-nodal-bnormals.plt").c_str(),std::ios::trunc);
  bool pointcloudout = false;
  f << "VARIABLES = x y " << (nbDim>2? 'z':' ')
               << " u v " << (nbDim>2? 'w':' ')
               << " a" << std::endl;
  for (CFuint i=0; i<btrs.size(); ++i) {

    const ConnectivityTable< CFuint > faces = *btrs[i]->getGeo2NodesConn();
    f << "ZONE DATAPACKING=POINT"
      << " ZONETYPE=" << (nbDim>2? "FETRIANGLE":"FELINESEG")
      << " N=" << nbNodes << " E=" << faces.nbRows()
      << (pointcloudout? (Ndim>2? " VARSHARELIST=([1-7]=1)" :
                                  " VARSHARELIST=([1-5]=1)") : "")
      << std::endl;
    if (!pointcloudout) {
      for (CFuint n=0; n<nbNodes; ++n)
        f << *h_nodes[n] << ' ' << h_mn_bnormal[n] << ' ' << h_mn_bnarea[n] << std::endl;
      pointcloudout = true;
    }
    for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc) {
      for (CFuint i=0; i<faces.nbCols(ifc); ++i)
        f << " " << faces(ifc,i)+1;
      f << std::endl;
    }

  }
  f.close();
  d.ver("output boundary normals.");
#endif
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setNodeVolume()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< CFreal > h_mn_volume = s_mn_volume.getDataHandle();

  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const CFuint nbNodes = h_nodes.size();
  const CFuint nbCells = geo2nodes->nbRows();

  // initialize
  h_mn_volume.resize(nbNodes);
  h_mn_volume = 0.;

  struct local_node_struct No_local[4];
  double vol;
  int inc_min;

  // elements volumes distributed to nodes
  for (CFuint ic=0; ic<nbCells; ++ic) {
    d.cellgeom(ic,No_local,&vol,&inc_min);
    for (int inc=0; inc<Nvtcell; ++inc)
      h_mn_volume[ No_local[inc].node ] += vol/(double) Nvtcell;
  }

  // get volume of partition
  double Vsub = 0.;
  for (CFuint n=0; n<nbNodes; ++n)
    if (h_nodes[n]->isParUpdatable())
      Vsub += h_mn_volume[n];

#if 0
  // synchronise nodal volumes
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  cf_assert(!d.m_restart);
  for (int n=0; n<nbNodes; ++n)
    (*h_states[n])[0] = h_mn_volume[n];
  d.synchronise();
  for (int n=0; n<nbNodes; ++n)
    h_mn_volume[n] = (*h_states[n])[0];
#endif

  // output partition volumes
  double Vtot = 0.;
  GlobalReduceOperation< GRO_SUM >(&Vsub,&Vtot);
  d.m_volume = Vtot;
  d.log("volume = " + StringOps::to_str((int) (100.*Vsub/Vtot)) + "% of " + StringOps::to_str(Vtot));


#if 0
  d.log("output nodal volume...");
  const CFuint nbDim = (h_nodes.size()?  h_nodes[0]->size() : 0);
  std::ofstream f(d.getFilename("debug-nodal-volume.plt").c_str(),std::ios::trunc);
  f << "VARIABLES = x y "<< (nbDim>2? 'z':' ') << " v" << std::endl;
  f << "ZONE DATAPACKING=POINT ZONETYPE=" << (nbDim>2? "FETETRAHEDRON":"FETRIANGLE") << " N=" << nbNodes << " E=" << nbCells << std::endl;
  for (CFuint n=0; n<nbNodes; ++n)
    f << *h_nodes[n] << ' ' << h_mn_volume[n] << std::endl;
  for (CFuint ic=0; ic<nbCells; ++ic) {
    for (int i=0; i<Nvtcell; ++i)
      f << " " << (*geo2nodes)(ic,i)+1;
    f << std::endl;
  }
  f.close();
  d.ver("output nodal volume.");
#endif
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setInitialSolution()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  // set VectorialFunction
  VectorialFunction f;
  f.setFunctions(d.m_initialvalues_def);
  f.setVariables(d.getNodalVariables());
  try {
    f.parse();
  }
  catch (ParserException& e) {
    d.log("VectorialFunction parsing: " + std::string(e.what()));
    throw;
  }

  // evaluate and set State's DataHandle
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  RealVector vars(0.,d.getNodalVariables().size());
  for (CFuint n=0; n<h_states.size(); ++n) {
    d.getNodalValues(n,vars);
    f.evaluate(vars,*h_states[n]);
  }

  // synchronise
  d.synchronise();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

