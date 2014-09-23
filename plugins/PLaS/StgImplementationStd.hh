#ifndef COOLFluiD_PLaS_StgImplementationStd_hh
#define COOLFluiD_PLaS_StgImplementationStd_hh

#include "Framework/MethodStrategy.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "PLaS/StgImplementation.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// This is a possible implementation of PLaS interface functions, specially
/// for simplex elements and Finite-Element type methods
class StgImplementationStd : public StgImplementation {

 public:

  /// Constructor
  StgImplementationStd(const std::string& name);

  /// Destructor
  ~StgImplementationStd() {}

  /// @returns the class name as a string
  static std::string getClassName() { return "StgImplementationStd"; }

  /// Sets up the implementation strategy
  void setup();


 private:  // calculate source sockets and boundaries

  /// Setup volumes statistics
  /// @param (return) number of elements in partition
  /// @param (return) partition volume
  /// @param (return) minimum element volume in partition
  /// @param (return) maximum element volume in partition
  void setupVolumes(int& Nelm, double& vtot, double& vmin, double& vmax);

  /// Setup mesh connectivities by setting pointer of element to nodes table
  /// @return pointer to inner cells connectivity table
  Common::SafePtr< Common::ConnectivityTable< CFuint > > setupConnectivity();

  /// Setup element to element connectivity table
  void setupElementToElement();


 private:  // sockets

  /// Socket to access nodes
  Framework::DataSocketSink< Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink< Framework::State*,Framework::GLOBAL > s_states;

  /// socket to access states (at the previous time-step)
  Framework::DataSocketSink< Framework::State* > s_paststates;

  /// Socket to access mapping surface face to corresponding boundary element
  Framework::DataSocketSink< std::pair< CFuint,CFuint > > s_faceneighcell;

  /// Socket to access nodal volumes
  Framework::DataSocketSink< CFreal > s_nvolume;

  /// Socket to access element volumes
  Framework::DataSocketSink< CFreal > s_evolume;

  /// Socket to access 'inner' elements normals
  Framework::DataSocketSink< std::vector< RealVector > > s_ienormals;

  /// Socket to access 'boundary' elements faces normals
  Framework::DataSocketSink< std::vector< RealVector > > s_benormals;


 public:  // sockets

  /// Returns DataSocket's this strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets() {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_nodes);
    r.push_back(&s_states);
    r.push_back(&s_paststates);
    r.push_back(&s_faceneighcell);
    r.push_back(&s_nvolume);
    r.push_back(&s_evolume);
    r.push_back(&s_ienormals);
    r.push_back(&s_benormals);
    return r;
  }


 public:  // data

  /// Number of dimensions (Node size)
  int nbDim;

  /// Connectivity list pointer, elements to nodes
  Common::SafePtr< Common::ConnectivityTable< CFuint > > m_geo2nodes;

  /// Connectivity list pointer, elements to elements
  Common::ConnectivityTable< int > m_geo2geo;


 public:  // PLaS interface, pure virtual functions implementations

  ////////////////////////////////////////////////////////////////////////////////////////
  //  Declaration of PLaS interface functions to be implemented in calling flow solver  //
  ////////////////////////////////////////////////////////////////////////////////////////

  /// Set flow solver parameters on initialization
  /// @param fp pointer to PLaS data structure
  void PLaS_SetFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);

  /// Set flow solver parameters on time step
  /// @param fp pointer to PLaS data structure
  void PLaS_SetFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp);

  /// Set partitioning data
  /// @param part partitioning data structure
  void PLaS_SetPartitioningData(PLAS_PART_DATA *part);

  /// Provide node of an element to PLaS
  /// @param elm element index
  /// @param enod node of the element
  /// @return node index
  int PLaS_GetElmNode(int elm, int enod);

  /// Provide neighbour of an element to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @return neighbour element index
  int PLaS_GetElmNeighbour(int elm, int eface);

  /// Provide coordinate of first node on a boundary face to PLaS
  /// @param bnd boundary index
  /// @param bface face of the boundary
  /// @param dim coordinate index
  /// @return coordinate of first node
  double PLaS_GetBndFaceRefCoord(int bnd, int bface, int dim);

  /// Provide domain element associated to a boundary face to PLaS
  /// @param bnd boundary index
  /// @param bface face of the boundary
  /// @param dummy not needed
  /// @return element index
  int PLaS_GetBndDomElm(int bnd, int bface, int dummy);

  /// Provide node coordinate to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return coordinate of the node
  double PLaS_GetNodCoord(int nod, int dim);

  /// Provide component of element normal to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @param dim coordinate index
  /// @return coordinate of the normal
  double PLaS_GetElmNormComp(int elm, int eface, int dim);

  /// Provide component of boundary face normal to PLaS
  /// @param elm element index
  /// @param bface face of the boundary
  /// @param dim coordinate index
  /// @return coordinate of the normal
  double PLaS_GetBndFaceNormComp(int bnd, int bface, int dim);

  /// Provide component of element face middle-point vector to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @param dim coordinate index
  /// @return coordinate of the face middle-point
  double PLaS_GetElmFaceMiddlePoint(int elm, int eface, int dim);

  /// Provide nodal area/volume to PLaS
  /// @param nod node index
  /// @return nodal cell volume
  double PLaS_GetNodVol(int nod);

  /// Provide element area/volume to PLaS
  /// @param elm element index
  /// @return element volume
  double PLaS_GetElmVol(int elm);

  /// Provide number of faces of a boundary to PLaS
  /// @param bnd boundary index
  /// @return number of boundary faces
  int PLaS_GetNumBndFaces(int bnd);

  /// Provide information about which boundary is a wall to PLaS
  /// @param bnd boundary index
  /// @return value if boundary is a wall, else zero
  int PLaS_GetWallBndFlag(int bnd);

  /// Provide information about which boundary is periodic to PLaS
  /// @param bnd boundary index
  /// @return value if boundary is periodic, else zero
  int PLaS_GetPerBndFlag(int bnd);

  /// Provide information about periodic boundary offset to PLaS
  /// @param bnd boundary index
  /// @param dim coordinate index
  /// @return periodic boundary offset
  double PLaS_GetPerBndOffset(int bnd, int dim);

  /// Provide type of an element to PLaS
  /// @param elm element index
  /// @return element type (codes in #defines)
  int PLaS_GetElementType(int elm);

  /// Provide nodal velocity component at time step n to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return velocity at time step n
  double PLaS_GetVelocityComp(int nod, int dim);

  /// Provide nodal velocity component at time step n-1 to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return velocity at time step n-1
  double PLaS_GetVelocityCompOld(int nod, int dim);

  /// Provide nodal temperature at time step n to PLaS
  /// @param nod node idx
  /// @return temperature at time step n
  double PLaS_GetTemperature(int nod);

  /// Provide nodal temperature at time step n-1 to PLaS
  /// @param nod = node idx
  /// @return temperature at time step n-1
  double PLaS_GetTemperatureOld(int nod);

  /// Provide nodal pressure at time step n to PLaS
  /// @param nod node idx
  /// @return pressure at time step n
  double PLaS_GetPressure(int nod);

  /// Provide nodal pressure at time step n-1 to PLaS
  /// @param nod node idx
  /// @return pressure at time step n-1
  double PLaS_GetPressureOld(int nod);

  /// Provide starting element for local brute force search
  /// @param pos position vector (xyz)
  /// @return lowest element for BF search (0)
  int PLaS_StartElementSearch(double *pos);

  /// Provide ending element for local brute force search
  /// @param pos position vector (xyz)
  /// @return highest element for BF search (numElm-1)
  int PLaS_EndElementSearch(double *pos);

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_PLaS_StgImplementationStd_hh

