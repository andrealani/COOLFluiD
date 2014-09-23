#ifndef COOLFluiD_Framework_ComputeFaceNormalsPyramP1_hh
#define COOLFluiD_Framework_ComputeFaceNormalsPyramP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "Framework/ComputeFaceNormalsFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

    /**
     * This class is a functor tha computes the face (outward) normals
     *
     * @author Andrea Lani
     */
class Framework_API ComputeFaceNormalsPyramP1 : public ComputeFaceNormalsFVMCC {
public:

  /**
   * Constructor
   */
  ComputeFaceNormalsPyramP1() : 
    ComputeFaceNormalsFVMCC(),
    _faceNormals(5,DIM_3D),
    _coordinates(5,DIM_3D)  
  {
  }
  
  /**
   * Destructor
   */
  ~ComputeFaceNormalsPyramP1()
  {
  }

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  void operator() (const CFuint& iFirstElem,
                   const CFuint& iLastElem,
                   const CFuint& iType);
  
private:
  
  RealMatrix _faceNormals;
  RealMatrix _coordinates;
  
}; // end of class ComputeFaceNormalsPyramP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeFaceNormalsPyramP1_hh
