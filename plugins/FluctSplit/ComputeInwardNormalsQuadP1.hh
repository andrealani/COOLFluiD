#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsQuadP1_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsQuadP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeInwardNormals.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor tha computes the face (outward) normals
/// @author Andrea Lani
/// @author Tiago Quintino

class FluctSplit_API ComputeInwardNormalsQuadP1 : public ComputeInwardNormals {

public:

  /// Constructor
  ComputeInwardNormalsQuadP1() : ComputeInwardNormals()
  {
  }

  /// Destructor
  ~ComputeInwardNormalsQuadP1()
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

}; // end of class ComputeInwardNormalsQuadP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsQuadP1_hh
