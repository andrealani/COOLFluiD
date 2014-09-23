#ifndef COOLFluiD_Framework_ComputeFaceNormalsQuadP1_hh
#define COOLFluiD_Framework_ComputeFaceNormalsQuadP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeFaceNormalsFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a functor that computes the face (outward) normals
 *
 * @author Andrea Lani
 *
 */
class Framework_API ComputeFaceNormalsQuadP1 : public ComputeFaceNormalsFVMCC {
public:

  /**
   * Constructor
   */
  ComputeFaceNormalsQuadP1() : ComputeFaceNormalsFVMCC()
  {
  }

  /**
   * Destructor
   */
  ~ComputeFaceNormalsQuadP1()
  {
  }

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  void operator() (const CFuint& iFirstElem,
		   const CFuint& iLastElem,
		   const CFuint& iType);

}; // end of class ComputeFaceNormalsQuadP1

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeFaceNormalsQuadP1_hh
