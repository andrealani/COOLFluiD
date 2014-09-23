#ifndef COOLFluiD_PLaS_StgImplementation_hh
#define COOLFluiD_PLaS_StgImplementation_hh

#include "PLaS/PLaSTrackingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// Strategy for implementation, naked pointer (defined in PLaSTrackingData)
extern StgImplementation* g_stg_implementation;

/**
 * This is a base implementation of PLaS interface functions, very empty
 * because it should only contain configuration and setup functions, leaving
 * the concrete PLaS interface implementation to derived classes
 */
class StgImplementation: public PLaSTrackingStg {

 public:

  /// This is the provider for all implementation strategies
  typedef Framework::BaseMethodStrategyProvider< PLaSTrackingData,StgImplementation > PROVIDER;

  /// Constructor
  StgImplementation(const std::string& name) : PLaSTrackingStg(name) {}

  /// Destructor
  virtual ~StgImplementation() {}

  /// Configures the object
  virtual void configure(Config::ConfigArgs& args) {
    PLaSTrackingStg::configure(args);
  }

  /// Sets up the object
  virtual void setup() {
    PLaSTrackingStg::setup();
  }

  /// Unsets the object
  virtual void unsetup() {
    PLaSTrackingStg::unsetup();
  }


 public:  // PLaS interface, pure virtual functions

  ////////////////////////////////////////////////////////////////////////////////////////
  //  Declaration of PLaS interface functions to be implemented in calling flow solver  //
  ////////////////////////////////////////////////////////////////////////////////////////

  /// Set flow solver parameters on initialization
  /// @param fp pointer to PLaS data structure
  virtual void PLaS_SetFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp) = 0;

  /// Set flow solver parameters on time step
  /// @param fp pointer to PLaS data structure
  virtual void PLaS_SetFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp) = 0;

  /// Set partitioning data
  /// @param part partitioning data structure
  virtual void PLaS_SetPartitioningData(PLAS_PART_DATA *part) = 0;

  /// Provide node of an element to PLaS
  /// @param elm element index
  /// @param enod node of the element
  /// @return node index
  virtual int PLaS_GetElmNode(int elm, int enod) = 0;

  /// Provide neighbour of an element to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @return neighbour element index
  virtual int PLaS_GetElmNeighbour(int elm, int eface) = 0;

  /// Provide coordinate of first node on a boundary face to PLaS
  /// @param bnd boundary index
  /// @param bface face of the boundary
  /// @param dim coordinate index
  /// @return coordinate of first node
  virtual double PLaS_GetBndFaceRefCoord(int bnd, int bface, int dim) = 0;

  /// Provide domain element associated to a boundary face to PLaS
  /// @param bnd boundary index
  /// @param bface face of the boundary
  /// @param dummy not needed
  /// @return element index
  virtual int PLaS_GetBndDomElm(int bnd, int bface, int dummy) = 0;

  /// Provide node coordinate to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return coordinate of the node
  virtual double PLaS_GetNodCoord(int nod, int dim) = 0;

  /// Provide component of element normal to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @param dim coordinate index
  /// @return coordinate of the normal
  virtual double PLaS_GetElmNormComp(int elm, int eface, int dim) = 0;

  /// Provide component of boundary face normal to PLaS
  /// @param elm element index
  /// @param bface face of the boundary
  /// @param dim coordinate index
  /// @return coordinate of the normal
  virtual double PLaS_GetBndFaceNormComp(int bnd, int bface, int dim) = 0;

  /// Provide component of element face middle-point vector to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @param dim coordinate index
  /// @return coordinate of the face middle-point
  virtual double PLaS_GetElmFaceMiddlePoint(int elm, int eface, int dim) = 0;

  /// Provide nodal area/volume to PLaS
  /// @param nod node index
  /// @return nodal cell volume
  virtual double PLaS_GetNodVol(int nod) = 0;

  /// Provide element area/volume to PLaS
  /// @param elm element index
  /// @return element volume
  virtual double PLaS_GetElmVol(int elm) = 0;

  /// Provide number of faces of a boundary to PLaS
  /// @param bnd boundary index
  /// @return number of boundary faces
  virtual int PLaS_GetNumBndFaces(int bnd) = 0;

  /// Provide information about which boundary is a wall to PLaS
  /// @param bnd boundary index
  /// @return value if boundary is a wall, else zero
  virtual int PLaS_GetWallBndFlag(int bnd) = 0;

  /// Provide information about which boundary is periodic to PLaS
  /// @param bnd boundary index
  /// @return value if boundary is periodic, else zero
  virtual int PLaS_GetPerBndFlag(int bnd) {
    return 0.;
  }

  /// Provide information about periodic boundary offset to PLaS
  /// @param bnd boundary index
  /// @param dim coordinate index
  /// @return periodic boundary offset
  virtual double PLaS_GetPerBndOffset(int bnd, int dim) {
    return 0.;
  }

  /// Provide type of an element to PLaS
  /// @param elm element index
  /// @return element type (codes in #defines)
  virtual int PLaS_GetElementType(int elm) = 0;

  /// Provide nodal velocity component at time step n to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return velocity at time step n
  virtual double PLaS_GetVelocityComp(int nod, int dim) {
    if (dim>=0 && dim<(int) getMethodData().m_vdef.size())
      return getMethodData().m_vdef[dim];
    return 0.;
  }

  /// Provide nodal velocity component at time step n-1 to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return velocity at time step n-1
  virtual double PLaS_GetVelocityCompOld(int nod, int dim) {
    return PLaS_GetVelocityComp(nod,dim);
  }

  /// Provide nodal velocity component derivative at time step n to PLaS
  /// @param nod node index
  /// @param idim variable to take derivatrive of
  /// @param jdim coordinate direction of derivative
  /// @return velocity derivative at time step n
  virtual double PLaS_GetVelocityDerivativeComp(int nod, int idim, int jdim) {
    return 0.;
  }

  /// Provide nodal velocity component derivative at time step n-1 to PLaS
  /// @param nod node index
  /// @param idim variable to take derivatrive of
  /// @param jdim coordinate direction of derivative
  /// @return velocity derivative at time step n-1
  virtual double PLaS_GetVelocityDerivativeCompOld(int nod, int idim, int jdim) {
    return PLaS_GetVelocityDerivativeComp(nod,idim,jdim);
  }

  /// Provide nodal temperature at time step n to PLaS
  /// @param nod node idx
  /// @return temperature at time step n
  virtual double PLaS_GetTemperature(int nod) {
    return getMethodData().m_tdef;
  }

  /// Provide nodal temperature at time step n-1 to PLaS
  /// @param nod = node idx
  /// @return temperature at time step n-1
  virtual double PLaS_GetTemperatureOld(int nod) {
    return PLaS_GetTemperature(nod);
  }

  /// Provide nodal pressure at time step n to PLaS
  /// @param nod node idx
  /// @return pressure at time step n
  virtual double PLaS_GetPressure(int nod) {
    return getMethodData().m_pdef;
  }

  /// Provide nodal pressure at time step n-1 to PLaS
  /// @param nod node idx
  /// @return pressure at time step n-1
  virtual double PLaS_GetPressureOld(int nod) {
    return PLaS_GetPressure(nod);
  }

  /// Provide starting element for local brute force search
  /// @param pos position vector (xyz)
  /// @return lowest element for BF search (0)
  virtual int PLaS_StartElementSearch(double *pos) = 0;

  /// Provide ending element for local brute force search
  /// @param pos position vector (xyz)
  /// @return highest element for BF search (numElm-1)
  virtual int PLaS_EndElementSearch(double *pos) = 0;

  /// Provide Eulerian time scale (tke/disspation) for a node
  /// @param nod node index
  /// @return Eulerian time scale (k/eps)
  virtual double PLaS_GetEulerianTimeScale(int nod) {
    return 1.;
  }


 public:  // PLaS interface, non-virtual functions

  /// Pass data for screen output
  /// @param text screen output
  void PLaS_ScreenOutput(char *text) { CFLog(INFO,"PLaS: " << text); }


  /// Pass data for a screen warning
  /// @param text screen warning
  void PLaS_ScreenWarning(char *text) { CFLog(WARN,"PLaS: " << text); }

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_PLaS_StgImplementation_hh

