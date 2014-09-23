#ifndef COOLFluiD_Numerics_FiniteVolume_KernelData_hh
#define COOLFluiD_Numerics_FiniteVolume_KernelData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores data to be exchanged with the GPU
 *
 * @author Andrea Lani
 *
 */
template <typename T> 
class KernelData {
public:
  
  /// Constructor
  HOST_DEVICE KernelData(CFuint nbCellsIn, 
			 T* statesIn, 
			 T* nodesIn,
			 T* centerNodesIn,
			 T* ghostStatesIn,
			 T* ghostNodesIn,
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
    nodes = nodesIn; 
    centerNodes = centerNodesIn;
    ghostStates = ghostStatesIn; 
    ghostNodes = ghostNodesIn; 
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
  T* nodes;
  T* centerNodes;
  T* ghostStates;
  T* ghostNodes;
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

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_KernelData_hh
