#ifndef COOLFluiD_Numerics_FluctSplit_InwardNormalsData_hh
#define COOLFluiD_Numerics_FluctSplit_InwardNormalsData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"

#include "Framework/PhysicalModel.hh"

#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class encapsulates data related with the normals to the faces,
/// allowing to map each normal to its coordinates and to the
/// corresponding face for FluctuationSplit method
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API InwardNormalsData {

public: // functions

  /// Constructor
  /// @param nodalNormals  the components of the all scaled inward normals
  /// @param faceNormals   the components of the all scaled face normals
  /// @param iType         type of the element
  /// @param nbFaces       number of faces
  /// @param nodalAreas    area (length in 2D) associated with a node
  /// @param faceAreas     area (length in 2D) associated with a face
  /// @pre                 all these data are computed in the concrete Connection,
  ///                      since this is the only object knowing what is required of
  ///                      the MeshData by the chosen SpaceMethod.
  InwardNormalsData(CFreal *const faceNormals,
  	    CFreal *const faceAreas,
  	    CFreal *const nodalNormals,
  	    CFreal *const nodalAreas,
  	    const CFuint nbFaces,
  	    const CFuint iType);
  /// Destructor
  ~InwardNormalsData();

  /// Get the pointer to the coordinates of the normal associated
  /// with a face.
  /// @param  iFace  local ID of the face in the Cell
  /// @return pointer to the first component of the normal
  const CFreal* getFaceNormalPtr(const CFuint iFace) const
  {
    cf_assert(iFace < m_nbFaces);
    return &m_faceNormals[iFace*m_dim];
  }

  /// Get the required coordinate of the normal associated
  /// with a dof
  /// @param  iState  local ID of the dof in the Cell
  /// @param  iCoord  ID of the space coordinate
  /// @return the specified space coordinate of the normal
  CFreal getNodalNormComp(const CFuint iNode,
                          const CFuint iCoord) const
  {
    if (isSimplex()) {
      return - m_scale * m_faceNormals[m_stateToFaceID[m_iType][iNode]*m_dim + iCoord];
    }
    else {
      return m_scale * m_nodalNormals[iNode*m_dim + iCoord];
    }
  }

  /// Get the required coordinate of the normal associated
  /// with a dof
  /// @param  iState  local ID of the dof in the Cell
  /// @param  iCoord  ID of the space coordinate
  /// @return the specified space coordinate of the normal
  CFreal getNodalUnitNormComp(const CFuint iNode,
                              const CFuint iCoord) const
  {
    return getNodalNormComp(iNode,iCoord) / getAreaNode(iNode);
  }

  /// Get the required coordinate of the normal associated
  /// with a face.
  /// @param  iFace  local ID of the face in the Cell
  /// @param  iCoord  ID of the space coordinate
  /// @return specified space coordinate of the normal
  CFreal getFaceNormComp(const CFuint iFace,
                         const CFuint iCoord) const
  {
    cf_assert(iFace < m_nbFaces);
    return m_scale * m_faceNormals[iFace*m_dim + iCoord];
  }

  /// Get the required coordinate of the normal associated
  /// with a face.
  /// @param  iFace  local ID of the face in the Cell
  /// @param  iCoord  ID of the space coordinate
  /// @return specified space coordinate of the normal
  CFreal getFaceUnitNormComp(const CFuint iFace,
                             const CFuint iCoord) const
  {
    return getFaceNormComp(iFace,iCoord) / getAreaFace(iFace);
  }

  /// Get the area of the specified face
  CFreal getAreaNode(const CFuint iNode) const
  {
    if (isSimplex()) {
      return std::abs(m_scale) * m_faceAreas[m_stateToFaceID[m_iType][iNode]];
    }
    else {
      return std::abs(m_scale) * m_nodalAreas[iNode];
    }
  }

  /// Get the area of the specified face
  CFreal getAreaFace(const CFuint iFace) const
  {
    cf_assert(iFace < m_nbFaces);
    return std::abs(m_scale) * m_faceAreas[iFace];
  }

  /// Get the type ID of the element
  CFuint getTypeID() const
  {
    return m_iType;
  }


  /// Set the normals of the specified face
  void setFaceNormals(const RealMatrix& p)
  {
    cf_assert(p.size() == m_nbFaces*m_dim);
    for (CFuint i = 0; i < p.size(); ++i) {
      m_faceNormals[i] = p[i];
    }
  }

  /// Set the area of the specified face
  void setNodalNormals(const RealMatrix& p)
  {
    for (CFuint i = 0; i < p.size(); ++i) {
      m_nodalNormals[i] = p[i];
    }
  }

  /// Set the area of the specified face
  void setNodalAreas(const RealVector& nodalAreas)
  {
    for (CFuint i = 0; i < nodalAreas.size(); ++i) {
      m_nodalAreas[i]  = nodalAreas[i];
    }
  }

  /// Set the area of the specified face
  void setFaceAreas(const RealVector& faceAreas)
  {
    cf_assert(faceAreas.size() == m_nbFaces);
    for (CFuint i = 0; i < faceAreas.size(); ++i) {
      m_faceAreas[i] = faceAreas[i];
    }
  }

  /// Set the type ID of the element
  void setTypeID(CFuint typeID)
  {
    m_iType = typeID;
  }

  /// @return number of faces
  CFuint nbFaces() const
  {
    return m_nbFaces;
  }

  /// Set the scale factor
  void scale(const CFreal& factor)
  {
    m_scale = factor;
  }

  /// Reset the scale factor to unity
  void unscale()
  {
    m_scale = 1.0;
  }

  /// Resize the array of state to opposite face IDs
  /// Also caches the dimension of the mesh
  /// @pre this has to be done in a place (Setup command) where the
  ///      number of element types is known
  static void resizeStateToFaceIDMap(const CFuint nbElemTypes)
  {
    m_stateToFaceID.resize(nbElemTypes);
    m_dim = Framework::PhysicalModelStack::getActive()->getDim();
  }

  /// Set the state to opposite face IDs for the given type
  static void setStateToFaceID(const CFuint iType,
  		       const std::vector<CFuint>& ids)
  {
    m_stateToFaceID[iType].resize(ids.size());
    m_stateToFaceID[iType] = ids;
  }

protected:

  bool isSimplex() const
  {
    return m_nodalNormals == CFNULL;
  }

private: // data

  /// ptrs to the face normals storage
  CFreal* m_faceNormals;

  /// array of the face areas
  CFreal* m_faceAreas;

  /// matrix storing all the nodal normals in the element
  CFreal* m_nodalNormals;

  /// list of the nodal areas
  CFreal* m_nodalAreas;

  /// number of faces
  CFuint  m_nbFaces;

  /// iType of this normal
  CFuint  m_iType;

  /// scale factor for normals
  CFreal  m_scale;

  /// mapping state to ID opposite face
  static std::vector< std::vector<CFuint> > m_stateToFaceID;

  /// cache of the number of dimensions of the mesh to avoid getting it all the time
  static CFuint m_dim;

}; // end of class InwardNormalsData

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_InwardNormalsData_hh
