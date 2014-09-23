#include "InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

// mapping state to ID opposite face
std::vector<std::vector<CFuint> > InwardNormalsData::m_stateToFaceID;

// cache for the dimension of the mesh
CFuint InwardNormalsData::m_dim = 0;

//////////////////////////////////////////////////////////////////////////////

InwardNormalsData::InwardNormalsData( CFreal *const faceNormals,
		                                  CFreal *const faceAreas,
		                                  CFreal *const nodalNormals,
		                                  CFreal *const nodalAreas,
		                                  const CFuint nbFaces,
		                                  const CFuint iType) :
  m_faceNormals(faceNormals),
  m_faceAreas(faceAreas),
  m_nodalNormals(nodalNormals),
  m_nodalAreas(nodalAreas),
  m_nbFaces(nbFaces),
  m_iType(iType),
  m_scale(1.0)
{
}

//////////////////////////////////////////////////////////////////////////////

InwardNormalsData::~InwardNormalsData()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
