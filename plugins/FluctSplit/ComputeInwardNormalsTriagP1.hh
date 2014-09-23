#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTriagP1_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTriagP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeInwardNormals.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor tha computes the face (outward) normals
/// @author Andrea Lani
/// @author Tiago Quintino

class FluctSplit_API ComputeInwardNormalsTriagP1 : public ComputeInwardNormals {
public:

  /// Constructor
  ComputeInwardNormalsTriagP1();

  /// Destructor
  ~ComputeInwardNormalsTriagP1()
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

}; // end of class ComputeInwardNormalsTriagP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsTriagP1_hh
