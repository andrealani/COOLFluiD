#ifndef COOLFluiD_Numerics_FiniteVolume_NullDerivativeComputer_hh
#define COOLFluiD_Numerics_FiniteVolume_NullDerivativeComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "DerivativeComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes derivatives using a diamond volume in 2D
 *
 * @author Andrea Lani
 *
 */
class NullDerivativeComputer : public DerivativeComputer {
public:

  /**
   * Constructor
   */
  NullDerivativeComputer(const std::string& name);

  /**
   * Default destructor
   */
  ~NullDerivativeComputer();

  /**
   * Set up the member data
   */
  virtual void setup();

  /*
   * Compute the gradients
   */
  void computeGradients(const RealMatrix& values,
			std::vector<RealVector*>& gradients);

  /**
   * Compute the control volume around the current face
   */
  void computeControlVolume(std::vector<RealVector*>& states, 
			    Framework::GeometricEntity *const geo);

  /**
   * Compute the average values corresponding to the given values
   */
  void computeAverageValues(Framework::GeometricEntity *const geo,
			    const std::vector<RealVector*>& values,
			    RealVector& avValues);
  
  /**
   * Get the maximum number of vertices in the control volume
   */
  CFuint getMaxNbVerticesInControlVolume() const
  {
    return 0;
  }

  /**
   * Get the current number of vertices in the control volume
   */
  CFuint getNbVerticesInControlVolume(Framework::GeometricEntity *const geo) const
  {
    return 0;
  }

  /**
   * Get the jacobian of the gradients
   */
  Common::SafePtr<std::vector<RealVector> > getGradientsJacob();

}; // end of class NullDerivativeComputer

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NullDerivativeComputer_hh
