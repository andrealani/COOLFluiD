#ifndef COOLFluiD_Numerics_FiniteVolume_StateDiffDerivative_hh
#define COOLFluiD_Numerics_FiniteVolume_StateDiffDerivative_hh

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
class StateDiffDerivative : public DerivativeComputer {
public:

  /**
   * Constructor
   */
  StateDiffDerivative(const std::string& name);

  /**
   * Default destructor
   */
  ~StateDiffDerivative();

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
  virtual void computeControlVolume(std::vector<RealVector*>& states, 
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
 
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private: // data
 
  /// storage of the cell volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// distance between the states LR
  CFreal _dr;

  /// direction vector for LR
  RealVector _vectorLR;

}; // end of class StateDiffDerivative

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StateDiffDerivative_hh
