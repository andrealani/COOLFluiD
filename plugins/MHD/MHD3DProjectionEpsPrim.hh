#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionEpsPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionEpsPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD3DProjectionEpsVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 3D for projection scheme
 * for primitive variables
 *
 * @author Andrea Lani
 */

class MHD3DProjectionEpsPrim : public MHD3DProjectionEpsVarSet {
public: //function

  /**
   * Constructor
   * @see MHD3DProjectionEps
   */
  MHD3DProjectionEpsPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~MHD3DProjectionEpsPrim();

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
   * Set the matrix of the right and left eigenvectors
   * and the matrix of the eigenvalues
   */
  void computeEigenValuesVectors(RealMatrix& rightEv,
                            RealMatrix& leftEv,
                            RealVector& eValues,
                            const RealVector& normal);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
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
   * @see EulerPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data,
			       Framework::State& state);


private: // data

}; // end of class MHD3DProjectionEpsPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Physics_MHD_MHD3DProjectionEpsPrim_hh
