#ifndef COOLFluiD_Numerics_FiniteVolume_ConstantPolyRecLin_hh
#define COOLFluiD_Numerics_FiniteVolume_ConstantPolyRecLin_hh

//////////////////////////////////////////////////////////////////////////////

#include "PolyReconstructorLin.hh"

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
class ConstantPolyRecLin : public PolyReconstructorLin {
public:

  /**
   * Constructor
   */
  ConstantPolyRecLin(const std::string& name);

  /**
   * Default destructor
   */
  ~ConstantPolyRecLin();

  /**
   * Compute the gredients
   */
  void computeGradients()
  {
    // no gradient is needed for first order reconstruction
  }

  /**
   * Compute the limiters for all the cells
   */
  void computeLimiters()
  {
    // no limiter is needed for first order reconstruction
  }

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
      PolyReconstructorLin::needsSockets();

    result.push_back(&socket_linearizedStates);
    result.push_back(&socket_linearizedGhostStates);

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

private: // member data

  /// socket for Linearized State's
  Framework::DataSocketSink<Framework::State*> socket_linearizedStates;

  /// socket for Linearized Ghost State's
  Framework::DataSocketSink<Framework::State*> socket_linearizedGhostStates;


}; // end of class ConstantPolyRecLin

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ConstantPolyRecLin_hh
