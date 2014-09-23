#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsPrismP1_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsPrismP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeInwardNormals.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor tha computes the face (outward) normals
/// @author Tiago Quintino
class FluctSplit_API ComputeInwardNormalsPrismP1 : public ComputeInwardNormals {

public:

  /// Constructor
  ComputeInwardNormalsPrismP1() : ComputeInwardNormals()
  {
  }

  /// Destructor
  ~ComputeInwardNormalsPrismP1()
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

}; // end of class ComputeInwardNormalsPrismP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormalsPrismP1_hh
