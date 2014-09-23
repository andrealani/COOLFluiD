#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTetraP2_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTetraP2_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeInwardNormals.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor tha computes the face (outward) normals
/// @author Andrea Lani
/// @author Tiago Quintino

class FluctSplit_API ComputeInwardNormalsTetraP2 : public ComputeInwardNormals {
public:

  /// Constructor
  ComputeInwardNormalsTetraP2();

  /// Destructor
  ~ComputeInwardNormalsTetraP2()
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

  /// Compute the averaged inward normals
  void average(const CFuint& iFirstCell,
               const CFuint& iLastCell,
               const CFuint& iType);

}; // end of class ComputeInwardNormalsTriagP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTriagP2_hh
