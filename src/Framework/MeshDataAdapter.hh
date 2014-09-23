// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MESHDATAADAPTER_HH
#define MESHDATAADAPTER_HH

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshData.hh"
#include "Framework/MeshDataInputSource.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class makes MeshData look like a MeshDataInputSource
/// It is primary goal is to provide non-performance critical
/// access to the raw mesh data. (Checkpointing, writing meshes, ...)
class Framework_API MeshDataAdapter : public MeshDataInputSource {
public:

  /// Constructor
  MeshDataAdapter () :
    _States(MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle()),
    _Nodes(MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle()),
    _Elements(MeshDataStack::getActive()->getTrs("InnerCells")) // rename "Elements" the old "InnerCells"
  {
  }

  /// Destructor
  virtual ~MeshDataAdapter ()
  {
  }

  virtual CFuint getDimension () const
  {
    return PhysicalModelStack::getActive()->getDim ();
  }

  virtual CFuint getNbEquations () const
  {
    return PhysicalModelStack::getActive()->getNbEq ();
  }

  virtual CFuint getNbNodes () const
  {
    return _Nodes.size ();
  }

  virtual CFuint getNbStates () const
  {
    return _States.size ();
  }


  virtual CFGeoEnt::Type getGeomType (CFuint CurTRS) const
  {
    cf_assert (CurTRS < getNbTRSs ());

    // geomtype looks unused...
    return CFGeoEnt::FACE;
  }

  virtual CFuint getNbElements () const
  {
    return _Elements->getLocalNbGeoEnts();
  }

  virtual CFuint getNbElementTypes () const
  {
    return MeshDataStack::getActive()->getElementTypeData()->size();
  }

  virtual CFPolyOrder::Type getGeometricPolyOrder () const
  {
    return static_cast<CFPolyOrder::Type>
      ((*MeshDataStack::getActive()->getElementTypeData())[0].getGeoOrder ());
  }

  virtual CFPolyOrder::Type getSolutionPolyOrder () const
  {
    return static_cast<CFPolyOrder::Type>
      ((*MeshDataStack::getActive()->getElementTypeData())[0].getSolOrder ());
  }

  /// maybe this function could go: we always return 0 anyway when we don't
  /// have a solution
  virtual bool isWithSolution () const
  {
    return true;
  }

  virtual CFuint getNbTRSs () const
  {
    return MeshDataStack::getActive()->getTrsList().size();
  }

  virtual CFuint getNbTRs (CFuint TRS) const
  {
    cf_assert (TRS < MeshDataAdapter::getNbTRSs ());
    return (MeshDataStack::getActive()->getTrsList()[TRS])->getNbTRs();
  }

  virtual std::string getNameTRS (CFuint TRS) const
  {
    cf_assert (TRS < MeshDataAdapter::getNbTRSs ());
    return MeshDataStack::getActive()->getTrsList()[TRS]->getName();
  }

  virtual CFuint getNbGeomEnts (CFuint TRS, CFuint TR) const
  {
    return (*MeshDataStack::getActive()->getTrsList()[TRS])[TR]->getLocalNbGeoEnts();
  }

  virtual void getElementType (CFuint Nb, ElementTypeData & Out) const
  {
    cf_assert (Nb < MeshDataAdapter::getNbElementTypes ());
    Out = (*MeshDataStack::getActive()->getElementTypeData())[Nb];
  }

  virtual void getGeoElement (Common::SafePtr<TopologicalRegionSet> trs,
                              CFuint TR,
                              CFuint Geom,
                              ElementDataType & Out) const
  {
    cf_assert (TR < trs->getNbTRs());
    cf_assert (Geom < (*trs)[TR]->getLocalNbGeoEnts());

    // this idx counts the geometric entities already processed
    // inside the current TRS
    const CFuint geoID = (*trs)[TR]->getGeoIDInTrs(Geom);
    Out.first.resize (trs->getNbNodesInGeo(geoID));
    Out.second.resize (trs->getNbStatesInGeo(geoID));

    for (CFuint i=0; i<Out.first.size(); ++i) {
      Out.first[i] = _Nodes[trs->getNodeID(geoID,i)]->getGlobalID();
    }

    // HACK for FVMCC
    const bool isFVMCC = (_Elements->getNbStatesInGeo(0) == 1);
    if (!isFVMCC) {
      for (CFuint i = 0; i < Out.second.size(); ++i) {
        Out.second[i] = _States[trs->getStateID(geoID,i)]->getGlobalID();
      }
    }
    else {
      Out.second.resize(1);
      Out.second[0] = _States[trs->getStateID(geoID,0)]->getGlobalID();
    }
  }

  virtual void getNode (CFuint NNode, NodeDataType & Out) const
  {
    cf_assert (NNode < _Nodes.size());
    const Node & N = *_Nodes[NNode];
    Out.resize(N.size());
    for (CFuint i=0; i<N.size(); ++i)
      Out[i]=N[i];
  }

  virtual CFuint getGlobalNodeID (CFuint LocalID) const
  {
    cf_assert (LocalID < _Nodes.size());
    return _Nodes[LocalID]->getGlobalID ();
  }

  virtual void getState (CFuint NState, StateDataType & Out) const
  {
    cf_assert (NState < _States.size());
    const State & S = *_States[NState];
    Out.resize (S.size());
    for (CFuint i=0; i<S.size(); ++i)
      Out[i]= S[i];
  }

  virtual CFuint getGlobalStateID (CFuint StateID) const
  {
    cf_assert (StateID < _States.size());
    return _States[StateID]->getGlobalID ();
  }

  virtual CFuint getGlobalElementID (CFuint ElementID) const
  {
    // shouldn't this be global ID ??
    cf_assert (ElementID < _Elements->getLocalNbGeoEnts());
    return _Elements->getGlobalGeoID(ElementID);
  }

  virtual void getElement (CFuint El, ElementDataType & Out) const
  {
    // shouldn't this be global ID ??
    cf_assert (El < _Elements->getLocalNbGeoEnts());

    const CFuint nbElemNodes = _Elements->getNbNodesInGeo(El);
    const CFuint nbElemStates = _Elements->getNbStatesInGeo(El);
    Out.first.resize (nbElemNodes);
    Out.second.resize (nbElemStates);

    // set the global IDs of the element nodes
    for (CFuint i = 0; i < nbElemNodes; ++i) {
      Out.first[i] = _Nodes[_Elements->getNodeID(El,i)]->getGlobalID();
    }

    // set the global IDs of the element states
    for (CFuint i = 0; i < nbElemStates; ++i) {
      Out.second[i] = _States[_Elements->getStateID(El,i)]->getGlobalID();
    }
  }

  std::vector< Common::SafePtr<TopologicalRegionSet> > getTrsList() const
  {
    return MeshDataStack::getActive()->getTrsList();
  }

  const std::vector<std::vector<CFuint> > & getTotalTRSInfo () const
  {
    return MeshDataStack::getActive()->getTotalTRSInfo ();
  }

  const std::vector<std::string> & getTotalTRSNames () const
  {
    return MeshDataStack::getActive()->getTotalTRSNames ();
  }

private:

  const DataHandle<State*,Framework::GLOBAL> _States;
  const DataHandle<Node*,Framework::GLOBAL>  _Nodes;
  const Common::SafePtr<TopologicalRegionSet> _Elements;

};

//////////////////////////////////////////////////////////////////////////////

    }
}

#endif
