#ifndef COOLFluiD_Physics_MHD_MHD2DCons_hh
#define COOLFluiD_Physics_MHD_MHD2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 */
class MHD2DCons : public MHD2DVarSet {
public: // classes

  /**
   * Constructor
   * @see MHD2D
   */
  MHD2DCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~MHD2DCons();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Get extra variable names
   */
  std::vector<std::string> getExtraVarNames() const;

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  void computeEigenValuesVectors(RealMatrix& rightEv,
                             RealMatrix& leftEv,
                             RealVector& eValues,
                             const RealVector& normal);

  /**
   * Set the rotation matrices
   */
  void setRotationMatrices(const RealVector& normals,
			   RealMatrix& rm,
			   RealMatrix& rmInv);

  /**
   * Get the state flux in rotated reference frame
   */
  RealVector& getRotatedFlux(const Framework::State& state,
		      const RealVector& normal);

  /**
   * Set the matrix of wave strengths
   */
  void setWaveStrengths(RealVector& waveStrengths,
			const RealMatrix& leftEv,
			Framework::State& leftState,
			Framework::State& rightState,
			const RealVector& normal);

  /**
   * Set the total magnetic field and energy values
   */
  void setDimensionalValuesPlusExtraValues(const Framework::State& state, 
					   RealVector& result,
					   RealVector& extra);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  void computePhysicalData(const Framework::State& state,
			   RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data,
			       Framework::State& state);
  
private: // data

  /// Rotation matrix
  RealMatrix _rm;

  /// Inverse of the rotation matrix
  RealMatrix _rmInv;

  /// temporary linear state
  RealVector _linearState;

  /// temporary linear state
  RealVector _linState;

  /// temporary state flux
  RealVector _stateFlux;

  /// Transformation matrix
  RealMatrix _dudw;

  /// Temporary matrix
  RealMatrix _rightMat;
  
  RealVector _leftStateVector;
  RealVector _rightStateVector;
  RealVector _leftStateVectorRot;
  RealVector _rightStateVectorRot;
  RealVector _dw;
  
}; // end of class MHD2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DCons_hh
