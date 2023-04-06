#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionEpsCons_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionEpsCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD3DProjectionEpsVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 3D for projection scheme
 * for conservative variables
 *
 * @author Andrea Lani
 */
class MHD3DProjectionEpsCons : public MHD3DProjectionEpsVarSet {
public: // classes

  /**
   * Constructor
   * @see MHD3DProjectionEps
   */
  MHD3DProjectionEpsCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~MHD3DProjectionEpsCons();

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
  RealVector& getRotatedFlux(const Framework::State& leftState,
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

}; // end of class MHD3DProjectionEpsCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionEpsCons_hh
