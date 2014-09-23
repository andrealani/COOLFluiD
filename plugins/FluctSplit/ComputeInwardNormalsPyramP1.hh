#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsPyramP1_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsPyramP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeInwardNormals.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor that computes the face (outward) normals
/// @author Tiago Quintino
class FluctSplit_API ComputeInwardNormalsPyramP1 : public ComputeInwardNormals {

public:

  /// Constructor
  ComputeInwardNormalsPyramP1() : ComputeInwardNormals()
  {
  }

  /// Destructor
  ~ComputeInwardNormalsPyramP1()
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


}; // end of class ComputeInwardNormalsPyramP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsPyramP1_hh
