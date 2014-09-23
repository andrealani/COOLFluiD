#ifndef COOLFluiD_Physics_MHD_MHD3DPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 3D for conservative
 * variables
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 */
class MHD3DPrim : public MHD3DVarSet {
public:

  /**
   * Constructor
   * @see MHD3D
   */
  MHD3DPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~MHD3DPrim();

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
   * Split the jacobian
   */
  void splitJacobian(RealMatrix& jacobPlus,
		  RealMatrix& jacobMin,
		  RealVector& eValues,
		  const RealVector& normal);

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  void computeEigenValuesVectors(RealMatrix& rightEv,
			     RealMatrix& leftEv,
			     RealVector& eValues,
			     const RealVector& normal);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see MHDPhysicalModel
   */
  void computePhysicalData(const Framework::State& state,
			   RealVector& data);

  /**
   * Set the total magnetic field and energy values
   */
  void setDimensionalValuesPlusExtraValues(const Framework::State& state,
                                           RealVector& result,
                                           RealVector& extra);

  /**
   * Set a State starting from the given PhysicalData
   * @see MHDPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data,
			       Framework::State& state);


private: // helper functions

  /**
   * Set the rotation matrices
   */
  void setRotationMatrices(const RealVector& normals);

private: // data
  /// Temporary matrix of right eigenvectors
  RealMatrix                          _tempR;

  /// Temporary matrix of left eigenvectors
  RealMatrix                          _tempL;

  /// Rotation matrix
  RealMatrix                          _rm;

  /// Rotation matrix
  RealMatrix                          _rmInv;

  /// temporary matrix of right eigenvalues
  RealMatrix                       _rightEv;

  /// temporary matrix of left eigenvalues
  RealMatrix                       _leftEv;

 }; // end of class MHD3DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DPrim_hh
