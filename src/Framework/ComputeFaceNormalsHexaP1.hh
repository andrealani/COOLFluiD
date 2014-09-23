#ifndef COOLFluiD_Framework_ComputeFaceNormalsHexaP1_hh
#define COOLFluiD_Framework_ComputeFaceNormalsHexaP1_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeFaceNormalsFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a functor tha computes the face (outward) normals
 *
 * @author Thomas Wuilbaut
 */
class Framework_API ComputeFaceNormalsHexaP1 : public ComputeFaceNormalsFVMCC {
public:

  /**
   * Constructor
   */
  ComputeFaceNormalsHexaP1() : ComputeFaceNormalsFVMCC()
  {
  }

  /**
   * Destructor
   */
  ~ComputeFaceNormalsHexaP1()
  {
  }

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  void operator() (const CFuint& iFirstElem,
                   const CFuint& iLastElem,
                   const CFuint& iType);

}; // end of class ComputeFaceNormalsHexaP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeFaceNormalsHexaP1_hh
