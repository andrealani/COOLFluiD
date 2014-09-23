#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsHexaP1_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsHexaP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeInwardNormals.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor tha computes the face (outward) normals
/// @author Tiago Quintino
class FluctSplit_API ComputeInwardNormalsHexaP1 : public ComputeInwardNormals {

public:

  /// Constructor
  ComputeInwardNormalsHexaP1() :
    ComputeInwardNormals(),
    _nodalNormals(8,DIM_3D),
    _nodalAreas(8),
    _faceNormals(6,DIM_3D),
    _faceAreas(6),
    _coordinates(8,DIM_3D),
    _tmp(DIM_3D)
  {
  }

  /// Destructor
  ~ComputeInwardNormalsHexaP1()
  {
  }

  /// Overloading of the operator () to make this class act as a
  /// functor
  void operator() (const CFuint& iFirstCell,
                   const CFuint& iLastCell,
                   const CFuint& iType);

  /// Update the inward normals
  void update(const CFuint& iFirstCell,
              const CFuint& iLastCell,
              const CFuint& iType);
private:

  RealMatrix _nodalNormals;

  RealVector _nodalAreas;

  RealMatrix _faceNormals;

  RealVector _faceAreas;

  RealMatrix _coordinates;

  RealVector _tmp;

}; // end of class ComputeInwardNormalsHexaP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsHexaP1_hh
