#ifndef COOLFluiD_Framework_ComputeFaceNormalsTriagP1_hh
#define COOLFluiD_Framework_ComputeFaceNormalsTriagP1_hh

//////////////////////////////////////////////////////////////////////////////

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
class Framework_API ComputeFaceNormalsTriagP1 : public ComputeFaceNormalsFVMCC {
  
public:
  
  /**
   * Constructor
   */
  ComputeFaceNormalsTriagP1() : ComputeFaceNormalsFVMCC()
  {
  }
  
  /**
   * Destructor
   */
  ~ComputeFaceNormalsTriagP1()
  {
  }
  
  /**
   * Overloading of the operator () to make this class act as a 
   * functor
   */
  void operator() (const CFuint& iFirstElem,
		   const CFuint& iLastElem,
		   const CFuint& iType);
  
}; // end of class ComputeFaceNormalsTriagP1
    
//////////////////////////////////////////////////////////////////////////////
      
  } // namespace Framework
  
} // namespace COOLFluiD 

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeFaceNormalsTriagP1_hh
