#ifndef COOLFluiD_FluxReconstructionMethod_KernelData_hh
#define COOLFluiD_FluxReconstructionMethod_KernelData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores data to be exchanged with the GPU for Flux Reconstruction
 *
 * @author Ray Vandenhoeck
 *
 */
template <typename T> 
class KernelData {
public:
  
  /// Constructor
  HOST_DEVICE KernelData(CFuint nbCellsIn, 
			 T* statesIn, 
			 T* updateCoeffIn, 
			 T* rhsIn,
			 T* normalsIn,
			 T* uxIn,
			 T* uyIn,
			 T* uzIn,
			 CFint* isOutwardIn)
  {
    nbCells = nbCellsIn; 
    states = statesIn; 
    updateCoeff = updateCoeffIn; 
    rhs = rhsIn;
    normals = normalsIn;
    uX = uxIn;
    uY = uyIn;
    uZ = uzIn;
    isOutward = isOutwardIn;
  }
  
  /// Destructor
  HOST_DEVICE ~KernelData() {}
  
  CFuint nbCells;
  T* states; 
  T* updateCoeff; 
  T* rhs;
  T* normals;
  T* uX;
  T* uY;
  T* uZ;
  CFint* isOutward;
};
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_KernelData_hh
