#ifndef COOLFluiD_Numerics_FiniteVolume_ConstantPolyRec_hh
#define COOLFluiD_Numerics_FiniteVolume_ConstantPolyRec_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_PolyRec.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a constant polynomial reconstructor for FVM
 *
 * @author Andrea Lani
 *
 */
class ConstantPolyRec : public FVMCC_PolyRec {
public:

  /**
   * Constructor
   */
  ConstantPolyRec(const std::string& name);

  /**
   * Default destructor
   */
  ~ConstantPolyRec();

  /**
   * Compute the gredients
   */
  void computeGradients()
  {
    // no gradient is needed for first order reconstruction
  }
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
      Framework::PolyReconstructor<CellCenterFVMData>::needsSockets();
    return result;
  }
  
private: // helper function
  
  /**
   * Extrapolate the solution in the face quadrature points
   */
  void extrapolateImpl(Framework::GeometricEntity* const face);

  /**
   * Extrapolate the solution in the face quadrature points
   */
  void extrapolateImpl(Framework::GeometricEntity* const face,
		       CFuint iVar, CFuint leftOrRight);
  
}; // end of class ConstantPolyRec

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ConstantPolyRec_hh
