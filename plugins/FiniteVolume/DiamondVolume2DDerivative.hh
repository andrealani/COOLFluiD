#ifndef COOLFluiD_Numerics_FiniteVolume_DiamondVolume2DDerivative_hh
#define COOLFluiD_Numerics_FiniteVolume_DiamondVolume2DDerivative_hh

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
class DiamondVolume2DDerivative : public DerivativeComputer {
public:

  /**
   * Constructor
   */
  DiamondVolume2DDerivative(const std::string& name);

  /**
   * Default destructor
   */
  ~DiamondVolume2DDerivative();

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
    return 4;
  }

  /**
   * Get the current number of vertices in the control volume
   */
  CFuint getNbVerticesInControlVolume(Framework::GeometricEntity *const geo) const
  {
    return 4;
  }

  /**
   * Get the jacobian of the gradients
   */
  Common::SafePtr<std::vector<RealVector> > getGradientsJacob();

private: // data

  /// normal 01 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n01;

  /// normal 12 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n12;

  /// normal 23 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n23;

  /// normal 30 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n30;

  /// quad control volume nodes
  std::vector<Framework::Node*> _nodes;

  /// mid node coordinates
  RealVector _midNode;

  /// array for the weights
  RealVector _weights;

}; // end of class DiamondVolume2DDerivative

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DiamondVolume2DDerivative_hh
