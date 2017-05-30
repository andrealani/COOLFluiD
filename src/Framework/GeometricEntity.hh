// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GeometricEntity_hh
#define COOLFluiD_Framework_GeometricEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include "Common/SharedPtr.hh"
#include "MathTools/RealMatrix.hh"

#include "Framework/State.hh"
#include "Framework/InterpolatorProperties.hh"
#include "Framework/CFGeoEnt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {
    class BaseGeometricEntityProvider;
    
//////////////////////////////////////////////////////////////////////////////

/// This class represents a GeometricEntity (Cell, Face, Edge, ...).
/// A GeometricEntity has a list of Node's (coordinates of the mesh
/// points) and a list of State's (state vectors in which the
/// solution is stored).
/// @see Node
/// @see State
/// It can also have neighbors, meaning with that a list of other geometric
/// entities connected in some way (known by the client code) to this
/// geometric entity. Having neighbors give the geometric entity the
/// possibility of having a knowledge of the "world" outside it.
/// As an example, the neighbors of a Cell can be Cell's or Face's
/// surrounding it.
/// @invariant _states.size > 0
/// @invariant _nodes.size > 0
/// @invariant _neighborGeos.size >= 0
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API GeometricEntity {
public:
  
  typedef BaseGeometricEntityProvider PROVIDER;
  
  /// Constructor
  /// @param states   list of the states in the geometric entity
  /// @param nodes  list of the nodes in the geometric entity
  /// @pre states.size() == 0
  /// @pre nodes.size() == 0
  /// @post _neighborGeos.size == 0.
  GeometricEntity();

  /// Pure virtual destructor (this makes this class not instantiatable)
  /// @pre all derived classes have to define a destructor
  virtual ~GeometricEntity() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "GeometricEntity";}
  
  /// Set the ID
  void setID(const CFuint geoID) {  _geoID = geoID; }

  /// Get the ID
  CFuint getID() const { return _geoID; }

  /// Resize Nodes
  void resizeNodes(CFuint nbNodes)
  {
    _nodes.resize(nbNodes);
  }

  /// Resize States
  void resizeStates(CFuint nbStates)
  {
    _states.resize(nbStates);
  }

  /// Set the State corresponding to the given local ID
  void setState(CFuint iState, State *const state)
  {
    cf_assert(iState < _states.size());
    _states[iState] = state;;
  }

  /// Set the Node corresponding to the given local ID
  void setNode(CFuint iNode, Node *const node)
  {
    cf_assert(iNode < _nodes.size());
    _nodes[iNode] = node;;
  }

  /// @return the list of states
  /// @post doesn't return NULL
  std::vector<State*>* getStates() { return &_states; }

  /// @return the list of states
  /// @post doesn't return NULL
  const std::vector<State*>& getStates() const { return _states; }

  /// @return the State corresponding to the given local
  ///         (in this GeometricEntity) ID
  /// @post doesn't return NULL
  State* getState(CFuint stateID) const  {  cf_assert(stateID < _states.size());  return _states[stateID]; }

  /// @return the Node corresponding to the given local
  ///         (in this GeometricEntity) ID
  /// @post doesn't return NULL
  Node* getNode(CFuint nodeID) const  {  cf_assert(nodeID < _nodes.size());  return _nodes[nodeID];  }

  /// Resize neighbor geos
  void resizeNeighborGeos(CFuint nbNeighbors) { _neighborGeos.resize(nbNeighbors); }

  /// @return the neighbor GeometricEntity corresponding to the
  /// given local (in this GeometricEntity) ID
  /// @post doesn't return NULL
  GeometricEntity* getNeighborGeo(CFuint iGeo) const
  {
    cf_assert(iGeo < _neighborGeos.size());
    return _neighborGeos[iGeo];
  }

  /// @return the list of nodes
  /// @post doesn't return NULL
  std::vector<Node*>* getNodes() {  return &_nodes; }

  /// @return the list of nodes
  /// @post doesn't return NULL
  const std::vector<Node*>& getNodes() const {  return _nodes; }

  /// @return the number of nodes
  CFuint nbNodes() const { return _nodes.size(); }

  /// @return the number of states
  CFuint nbStates() const  {  return _states.size(); }

  /// @return the number of neighbouring geos
  CFuint nbNeighborGeos() const
  {
    return _neighborGeos.size();
  }

  /// Reset the list of State's
  void resetStates(const std::vector<State*>& states)
  {
    _states = states;
  }

  /// Reset the list of Node's
  void resetNodes(const std::vector<Node*>& nodes)
  {
    _nodes = nodes;
  }

  /// Add a State to the list of State's
  void addState(State* const state)
  {
    _states.push_back(state);
  }

  /// Add a Node to the list of Node's
  void addNode(Node* const node)
  {
    _nodes.push_back(node);
  }

    /// Set the neighbors of the current geometric entity.
    /// @param  neighborGeos list of the neighbor geometric entity
    ///         of the current one
  void setNeighborGeos(const std::vector<GeometricEntity*>& neighborGeos)
  {
    _neighborGeos = neighborGeos;
  }

    /// Set one neighbor geometric entity of the current the geometric entity.
    /// @param  neighborGeo  neighbor geometric entity
  void setNeighborGeo(CFuint iGeo, GeometricEntity *const neighborGeo)
  {
    cf_assert(iGeo < _neighborGeos.size());
    _neighborGeos[iGeo] = neighborGeo;
  }

  /// Say if the neighbor geometrical entities have already been set
  /// @pre _neighborGeos.size() == 0 in the constructor
  /// @return true if set, false if not
  bool areNeighborGeosSet() const
  {
    return (_neighborGeos.size() > 0) ? true : false;
  }

  /// Computes the volume, area or length
  /// @pre Are considered volume of a geometric entity:
  ///      in 3D: cell volume, face area, edge length
  ///      in 2D: cell area, face length
  ///      in 1D: cell length
  /// @todo will no longer b pure virtual
  ///        getVolumeIntegrator().integrateConstantFunctorOnGeo<ConstantFunctor<DIM_1D> >(cell,functor,result);
  virtual CFreal computeVolume() = 0;


  /// Computes the centroid of the geometric entity
  virtual RealVector computeCentroid() = 0;

  /// Computes maped coordinates of states
  virtual void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords) = 0;

  /// Get the neighbor geometric entities
  /// @pre Are considered neighbor geometric entities:
  ///      for a Cell:  the list of its Faces
  ///      for a Face:  the two Cells sharing that Face
  ///      for an Edge: the list of Faces sharing that edge
  const std::vector<GeometricEntity*>* getNeighborGeos() const
  {
    return &_neighborGeos;
  }

  /// Check if the given list of nodes matches the list of
  /// nodes of the geometric entity and counts the nodes that
  /// matches
  /// @return the number of matching nodes
  inline CFuint countMatchingNodes(const std::vector<Node*>& nodes);

  /// Check if the geometric entity contains a given state
  /// @return true or false
  inline bool containState(const State* state);

  /// Check if the geometric entity contains a given state
  /// @post if the state is found, its local ID will be set in iState
  /// @return true if the given state is found
  inline bool containState(const State* state, CFuint& iState);

  /// Check if the geometric entity contains a given geometric entity
  /// as neighbor
  /// @post if the geo is found, its local ID will be set in iGeo
  /// @return true if the given geo is found
  inline bool containGeo(GeometricEntity* geo, CFuint& iGeo);

  /// Check if the geometric entity contains a given node
  /// @return true or false
  inline bool containNode(const Node* node);

  /// Check if the geometric entity contains a given node
  /// @post if the node is found, its local ID will be set in iNode
  /// @return true if the given node is found
  inline bool containNode(const Node* node, CFuint& iNode);
  
  /// Checks if a given GeometricEntity shares states with this one
  inline bool sharesStatesWith(GeometricEntity& other);
  
  /// Checks if a given GeometricEntity shares nodes with this one
  inline bool sharesNodesWith(GeometricEntity& other);
  
  /// Computes the gradients of the shapefunctions at given coordinates
  virtual std::vector<RealMatrix> computeSolutionShapeFunctionGradients(const std::vector<RealVector>& coord) = 0;

  /// Computes the gradients of the shapefunctions at given coordinates
  virtual std::vector<RealMatrix> computeGeoShapeFunctionGradients(const std::vector<RealVector>& coord) = 0;

  /// Computes the gradients of the shapefunctions at given mapped coordinates
  virtual std::vector<RealMatrix> computeSolutionShapeFunctionGradientsInMappedCoordinates(const std::vector<RealVector>& mappedCoord) const = 0;

  /// Computes the Jacobian matrix of the geometric shapefunctions at given mapped coordinates
  virtual std::vector<RealMatrix> computeGeometricShapeFunctionJacobianMatrix(const std::vector<RealVector>& mappedCoord) const = 0;

  /// Computes the Jacobian determinant of the geometric shapefunctions at given mapped coordinates
  virtual std::valarray<CFreal> computeGeometricShapeFunctionJacobianDeterminant(const std::vector<RealVector>& mappedCoord) const = 0;

  /// Computes the value of the shapefunctions at given coordinates
  virtual RealVector computeShapeFunctionAtCoord(const RealVector& coord) = 0;

  /// Computes the value of the shapefunctions at given coordinates
  virtual RealVector computeGeoShapeFunctionAtCoord(const RealVector& coord) = 0;

  /// Computes the coordinates from the mapped coordinates
  virtual RealVector computeCoordFromMappedCoord(const RealVector& coord) = 0;

  /// Computes the mapped coordinates from the coordinates
  virtual RealVector computeMappedCoordFromCoord(const RealVector& coord) = 0;

  /// Computes the geometric shape functions from the mapped coordinates
  virtual RealVector computeGeoShapeFunctionAtMappedCoord(const RealVector& mappedCoord) = 0;

  /// Computes the shape functions from the mapped coordinates
  virtual RealVector computeShapeFunctionAtMappedCoord(const RealVector& mappedCoord) = 0;

  /// Computes the average normals of the face
  /// (Normals have the same dimensionality as the Cell)
  virtual std::vector<RealVector> computeAvgFaceNormals() = 0;

  /// Computes the normal of the face at given coordinates
  /// (Normals have the same dimensionality as the Cell)
  virtual std::vector<RealVector> computeFaceNormals(const RealVector& coord) = 0;

  /// Computes the normal to a given mapped coordinate plane, at the given mapped coordinates,
  /// scaled with the 'face Jacobian determinant'.
  /// (Normal has the same dimensionality as the Cell)
  virtual std::vector< RealVector >
      computeMappedCoordPlaneNormalAtMappedCoords(const std::vector< CFuint >& planeIdx,
                                                  const std::vector< RealVector >& coord) = 0;

  /// Computes the normal to a face at the given mapped coordinates,
  /// scaled with the 'face Jacobian determinant'.
  /// (Normal has the dimensionality of the Face + 1)
  virtual std::vector< RealVector >
      computeFaceJacobDetVectorAtMappedCoords(const std::vector< RealVector >& coord) = 0;

  /// Computes the average normal of the Cell
  /// (Normal has the dimensionality of the Cell + 1)
  virtual RealVector computeAvgCellNormal() = 0;

  /// Computes the normal of the Cell at given coordinates
  /// (Normal has the dimensionality of the Cell + 1)
  virtual RealVector computeCellNormal(const RealVector& coord) = 0;

  /// Gets the number of nodes of the solution shape function
  virtual CFuint getNbNodesSolutionShapeFunction() const = 0;

  /// Gets the number of nodes of the geometry shape function
  virtual CFuint getNbNodesGeometryShapeFunction() const = 0;

  /// Gets the ID for the solution interpolation
  virtual InterpolatorID getSolInterpolatorID() const = 0;

  /// Gets the ID for the geometric interpolation
  virtual InterpolatorID getGeomInterpolatorID() const = 0;

  /// Get the number of facets
  virtual CFuint getNbFacets() const = 0;

  /// Get the inheritant dimensionality of the GeometricEntity
  virtual CFuint getDimensionality() const = 0;

  /// Gets the CFGeoEnt::Type of the Geometric Entity
  virtual CFGeoEnt::Type getGeometryType() const = 0;

  /// Gets the type of CFGeoShape::Type
  virtual CFGeoShape::Type getShape() const = 0;

  /// Gets the type of geometry shape function
  virtual std::string getGeometryShapeFunctionName() const = 0;

  /// Gets the type of geometry shape function
  virtual CFPolyForm::Type getGeometryShapeFunctionType() const = 0;

  /// Gets the solution shape function order
  virtual CFPolyOrder::Type getGeometryShapeFunctionOrder() const = 0;

  /// Gets the type of geometry shape function
  virtual std::string getSolutionShapeFunctionName() const = 0;

  /// Gets the type of solution shape function
  virtual CFPolyForm::Type getSolutionShapeFunctionType() const = 0;

  /// Gets the solution shape function  order
  virtual CFPolyOrder::Type getSolutionShapeFunctionOrder() const = 0;

  /// Check if a point (defined with mapped coordinates) is inside an element
  virtual bool isInMappedElement(const RealVector& mappedCoord) = 0;

  /// Check if a point is inside an element
  virtual bool isInElement(const RealVector& coord) = 0;

  /// Release this GeometricEntity allowing to create it again
  void release() { _isUsed = false; }

protected: //data

  /// global ID of this GeometricEntity
  /// cannot be treated via IndexList because if you have
  /// 2 different sets of GeometricEntity, each one must have an
  /// independent numbering
  CFuint _geoID;

  /// flag telling if this GeometricEntity is used
  bool _isUsed;

  /// acquaintance of the states present in this GeometricEntity
  std::vector<State*> _states;

  /// acquaintance of the nodes present in this GeometricEntity
  std::vector<Node*> _nodes;

  /// neighbor GeometricEntity's : they can be different concrete
  /// entities (Cell's, Face's, Edge's) depending on the numerical
  /// method that is used
  std::vector<GeometricEntity*> _neighborGeos;

}; // end of class GeometricEntity

//////////////////////////////////////////////////////////////////////////////

/// definition of a GeometricEntity List
typedef std::vector<GeometricEntity*> GeomEntList;

/// definition of a GeometricEntity pointer with ownership
typedef Common::SharedPtr<GeometricEntity> GeomEntityPtr;

//////////////////////////////////////////////////////////////////////////////

inline CFuint GeometricEntity::
  countMatchingNodes(const std::vector<Node*>& nodes)
{
  CFuint matchingNodes = 0;
  const CFuint nbNodesPerFace = nodes.size();

  // loop over the nodes of the face
  for (CFuint ik = 0; ik < nbNodesPerFace; ++ik) {

    // loop over the nodes provided
    for (CFuint in = 0; in < _nodes.size(); ++in) {

      // add to the counter if match
      // and break out of inner loop
      if (_nodes[in] == nodes[ik]) {
        ++matchingNodes;
        break;
      }
    }
  }

  return matchingNodes;
}

//////////////////////////////////////////////////////////////////////////////

inline bool GeometricEntity::containState(const State* state)
{
  std::vector<State*>::iterator it;
  for (it = _states.begin(); it != _states.end(); ++it) {
    if(*it == state) return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

inline bool GeometricEntity::containState(const State* state,
                                            CFuint& iState)
{
  iState = 0;
  std::vector<State*>::iterator it;
  for (it = _states.begin(); it != _states.end(); ++it) {
    if (*it == state) return true;
    ++iState;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

inline bool GeometricEntity::containGeo(GeometricEntity* geo,
                              CFuint& iGeo)
{
  iGeo = 0;
  std::vector<GeometricEntity*>::iterator it;
  for (it = _neighborGeos.begin(); it != _neighborGeos.end(); ++it) {
    if (*it == geo) return true;
    ++iGeo;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

inline bool GeometricEntity::containNode(const Node* node)
{
  std::vector<Node*>::iterator it;
  for (it = _nodes.begin(); it != _nodes.end(); ++it) {
    if(*it == node) return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

inline bool GeometricEntity::containNode(const Node* node, CFuint& iNode)
{
  iNode = 0;
  std::vector<Node*>::iterator it;
  for (it = _nodes.begin(); it != _nodes.end(); ++it) {
    if (*it == node) return true;
    ++iNode;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool GeometricEntity::sharesStatesWith(GeometricEntity& other)
{
  using namespace std;

  const vector<State*>* states = other.getStates();
  cf_assert(states != CFNULL);

  CFuint nbStates = states->size();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    if(containState((*states)[iState])) return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool GeometricEntity::sharesNodesWith(GeometricEntity& other)
{
  using namespace std;
  
  const vector<Node*>* nodes = other.getNodes();
  cf_assert(nodes != CFNULL);

  CFuint nbNodes = nodes->size();
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    if(containNode((*nodes)[iNode])) return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GeometricEntity_hh
