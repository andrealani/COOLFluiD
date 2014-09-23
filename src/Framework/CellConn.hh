#ifndef COOLFluiD_Framework_CellConn_hh
#define COOLFluiD_Framework_CellConn_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
      
//////////////////////////////////////////////////////////////////////////////

class Framework_API CellConn {
public:
  
  /// set the number of faces in the cell
  /// @param nbFaces number of faces
  HOST_DEVICE void setNbFaces(CFuint nbFaces) {m_nbFaces = nbFaces;}
  
  /// get the number of faces in the cell
  HOST_DEVICE CFuint getNbFaces() const {return m_nbFaces;}
  
  /// set the number of nodes for the given face
  /// @param iFace        local (inside the cell) face ID
  /// @param nbFaceNodes  number of nodes in the given face
  HOST_DEVICE void setNbFaceNodes(CFuint iFace, CFuint nbFaceNodes) 
  {m_nbFaceNodes[iFace] = nbFaceNodes;}
  
  /// get the number of nodes for the given face
  /// @param iFace  local (inside the cell) face ID
  HOST_DEVICE CFuint getNbFaceNodes(CFuint iFace) const {return m_nbFaceNodes[iFace];}
  
  /// set the requested node in the given face
  /// @param iFace   local (inside the cell) face ID
  /// @param iNode   local (inside the face) node ID
  /// @param nodeID  local (inside the cell) node ID
  HOST_DEVICE void setNodeID(CFuint iFace, CFuint iNode, CFuint nodeID) {m_faceNodeConn[iFace][iNode] = nodeID;}
  
  /// get the requested node in the given face
  /// @param iFace  local (inside the cell) face ID
  /// @param iNode  local (inside the face) node ID
  HOST_DEVICE CFuint getNodeID(CFuint iFace, CFuint iNode) const {return m_faceNodeConn[iFace][iNode];}
  
private:
  
  /// number of faces
  CFuint m_nbFaces;
  
  /// number of nodes in faces (over-dimensioned using hexahedra as limiting example)
  CFuint m_nbFaceNodes[6];
  
  /// face-node connectivity (over-dimensioned using hexahedra as limiting example)
  CFuint m_faceNodeConn[6][4];
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CellConn_hh
