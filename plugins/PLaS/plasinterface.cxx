
/*
 * Implementation of PLaS function calls: it's a layer to convert the C calls
 * to the StgImplementation strategy calls of similar name. Everytime the
 * library interface changes this file has to be updated and so has the "real"
 * implementation file (the strategy) as well.
 */


#include "PLaS/StgImplementation.hh"


namespace COOLFluiD {
  namespace PLaS {


// PLaS functions delegation to implementation strategy
extern "C" {

  /// Set flow solver parameters on initialization
  /// @param fp pointer to PLaS data structure
  void plasinterface_setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp) {
    g_stg_implementation->PLaS_SetFlowSolverParamOnInit(fp);
  }


  /// Set flow solver parameters on time step
  /// @param fp pointer to PLaS data structure
  void plasinterface_setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp) {
    g_stg_implementation->PLaS_SetFlowSolverParamOnTimeStep(fp);
  }


  /// Set partitioning data
  /// @param part partitioning data structure
  void plasinterface_setPartitioningData(PLAS_PART_DATA *part) {
    g_stg_implementation->PLaS_SetPartitioningData(part);
  }


  /// Provide node of an element to PLaS
  /// @param elm element index
  /// @param enod node of the element
  /// @return node index
  int plasinterface_getElmNode(int elm, int enod) {
    return g_stg_implementation->PLaS_GetElmNode(elm,enod);
  }


  /// Provide neighbour of an element to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @return neighbour element index
  int plasinterface_getElmNeighbour(int elm, int eface) {
    return g_stg_implementation->PLaS_GetElmNeighbour(elm,eface);
  }


  /// Provide coordinate of first node on a boundary face to PLaS
  /// @param bnd boundary index
  /// @param bface face of the boundary
  /// @param dim coordinate index
  /// @return coordinate of first node
  double plasinterface_getBndFaceRefCoord(int bnd, int bface, int dim) {
    return g_stg_implementation->PLaS_GetBndFaceRefCoord(bnd,bface,dim);
  }


  /// Provide domain element associated to a boundary face to PLaS
  /// @param bnd boundary index
  /// @param bface face of the boundary
  /// @param dummy not needed
  /// @return element index
  int plasinterface_getBndDomElm(int bnd, int bface, int dummy) {
    return g_stg_implementation->PLaS_GetBndDomElm(bnd,bface,dummy);
  }


  /// Provide node coordinate to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return coordinate of the node
  double plasinterface_getNodCoord(int nod, int dim) {
    return g_stg_implementation->PLaS_GetNodCoord(nod,dim);
  }


  /// Provide component of element normal to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @param dim coordinate index
  /// @return coordinate of the normal
  double plasinterface_getElmNormComp(int elm, int eface, int dim) {
    return g_stg_implementation->PLaS_GetElmNormComp(elm,eface,dim);
  }


  /// Provide component of boundary face normal to PLaS
  /// @param elm element index
  /// @param bface face of the boundary
  /// @param dim coordinate index
  /// @return coordinate of the normal
  double plasinterface_getBndFaceNormComp(int bnd, int bface, int dim) {
    return g_stg_implementation->PLaS_GetBndFaceNormComp(bnd,bface,dim);
  }


  /// Provide component of element face middle-point vector to PLaS
  /// @param elm element index
  /// @param eface face of the element
  /// @param dim coordinate index
  /// @return coordinate of the face middle-point
  double plasinterface_getElmFaceMiddlePoint(int elm, int eface, int dim) {
    return g_stg_implementation->PLaS_GetElmFaceMiddlePoint(elm,eface,dim);
  }


  /// Provide nodal area/volume to PLaS
  /// @param nod node index
  /// @return nodal cell volume
  double plasinterface_getNodVol(int nod) {
    return g_stg_implementation->PLaS_GetNodVol(nod);
  }


  /// Provide element area/volume to PLaS
  /// @param elm element index
  /// @return element volume
  double plasinterface_getElmVol(int elm) {
    return g_stg_implementation->PLaS_GetElmVol(elm);
  }


  /// Provide number of faces of a boundary to PLaS
  /// @param bnd boundary index
  /// @return number of boundary faces
  int plasinterface_getNumBndFaces(int bnd) {
    return g_stg_implementation->PLaS_GetNumBndFaces(bnd);
  }


  /// Provide information about which boundary is a wall to PLaS
  /// @param bnd boundary index
  /// @return value if boundary is a wall, else zero
  int plasinterface_getWallBndFlag(int bnd) {
    return g_stg_implementation->PLaS_GetWallBndFlag(bnd);
  }


  /// Provide information about which boundary is periodic to PLaS
  /// @param bnd boundary index
  /// @return value if boundary is periodic, else zero
  int plasinterface_getPerBndFlag(int bnd) {
    return g_stg_implementation->PLaS_GetPerBndFlag(bnd);
  }


  /// Provide information about periodic boundary offset to PLaS
  /// @param bnd boundary index
  /// @param dim coordinate index
  /// @return periodic boundary offset
  double plasinterface_getPerBndOffset(int bnd, int dim) {
    return g_stg_implementation->PLaS_GetPerBndOffset(bnd,dim);
  }


  /// Provide type of an element to PLaS
  /// @param elm element index
  /// @return element type (codes in #defines)
  int plasinterface_getElementType(int elm) {
    return g_stg_implementation->PLaS_GetElementType(elm);
  }


  /// Provide nodal velocity component at time step n to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return velocity at time step n
  double plasinterface_getVelocityComp(int nod, int dim) {
    return g_stg_implementation->PLaS_GetVelocityComp(nod,dim);
  }


  /// Provide nodal velocity component at time step n-1 to PLaS
  /// @param nod node index
  /// @param dim coordinate index
  /// @return velocity at time step n-1
  double plasinterface_getVelocityCompOld(int nod, int dim) {
    return g_stg_implementation->PLaS_GetVelocityCompOld(nod,dim);
  }


  /// Provide nodal velocity component derivative at time step n to PLaS
  /// @param nod node index
  /// @param idim variable to take derivatrive of
  /// @param jdim coordinate direction of derivative
  /// @return velocity derivative at time step n-1
  double plasinterface_getVelocityDerivativeComp(int nod, int idim, int jdim) {
    return g_stg_implementation->PLaS_GetVelocityDerivativeComp(nod,idim,jdim);
  }


  /// Provide nodal velocity component derivative at time step n-1 to PLaS
  /// @param nod node index
  /// @param idim variable to take derivatrive of
  /// @param jdim coordinate direction of derivative
  /// @return velocity derivative at time step n-1
  double plasinterface_getVelocityDerivativeCompOld(int nod, int idim, int jdim) {
    return g_stg_implementation->PLaS_GetVelocityDerivativeCompOld(nod,idim,jdim);
  }


  /// Provide nodal temperature at time step n to PLaS
  /// @param nod node idx
  /// @return temperature at time step n
  double plasinterface_getTemperature(int nod) {
    return g_stg_implementation->PLaS_GetTemperature(nod);
  }


  /// Provide nodal temperature at time step n-1 to PLaS
  /// @param nod = node idx
  /// @return temperature at time step n-1
  double plasinterface_getTemperatureOld(int nod) {
    return g_stg_implementation->PLaS_GetTemperatureOld(nod);
  }


  /// Provide nodal pressure at time step n to PLaS
  /// @param nod node idx
  /// @return pressure at time step n
  double plasinterface_getPressure(int nod) {
    return g_stg_implementation->PLaS_GetPressure(nod);
  }


  /// Provide nodal pressure at time step n-1 to PLaS
  /// @param nod node idx
  /// @return pressure at time step n-1
  double plasinterface_getPressureOld(int nod) {
    return g_stg_implementation->PLaS_GetPressureOld(nod);
  }


  /// Provide starting element for local brute force search
  /// @param pos position vector (xyz)
  /// @return lowest element for BF search (0)
  int plasinterface_StartElementSearch(double *pos) {
    return g_stg_implementation->PLaS_StartElementSearch(pos);
  }


  /// Provide ending element for local brute force search
  /// @param pos position vector (xyz)
  /// @return highest element for BF search (numElm-1)
  int plasinterface_EndElementSearch(double *pos) {
    return g_stg_implementation->PLaS_EndElementSearch(pos);
  }


  /// Provide Eulerian time scale (tke/disspation) for a node
  /// @param nod node index
  /// @return Eulerian time scale (k/eps)
  double plasinterface_getEulerianTimeScale(int nod) {
    return g_stg_implementation->PLaS_GetEulerianTimeScale(nod);
  }


  /// Pass data for screen output
  /// @param text screen output
  void plasinterface_screenOutput(char *text) {
    g_stg_implementation->PLaS_ScreenOutput(text);
  }


  /// Pass data for a screen warning
  /// @param text screen warning
  void plasinterface_screenWarning(char *text) {
    g_stg_implementation->PLaS_ScreenWarning(text);
  }

}


  }  // namespace PLaS
}  // namespace COOLFluiD

