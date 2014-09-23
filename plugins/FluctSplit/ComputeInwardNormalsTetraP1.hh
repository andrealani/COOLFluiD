#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTetraP1_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTetraP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeInwardNormals.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor tha computes the face (outward) normals
/// @author Andrea Lani

class FluctSplit_API ComputeInwardNormalsTetraP1 : public ComputeInwardNormals {
public:

  /// Constructor
  ComputeInwardNormalsTetraP1();

  /// Destructor
  ~ComputeInwardNormalsTetraP1()
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

  /// average the inward normals
  void average(const CFuint& iFirstCell,
               const CFuint& iLastCell,
               const CFuint& iType);

}; // end of class ComputeInwardNormalsTetraP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTetraP1_hh
