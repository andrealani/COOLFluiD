#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionPrimE_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionPrimE_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD3DProjectionPrim.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 3D for projection scheme
 * for primitive variables with decomposition of the energy equation
 *
 * @author Andrea Lani
 * @author Haopeng Wang
 */
class MHD3DProjectionPrimE : public MHD3DProjectionPrim {
public: //function

  /**
   * Constructor
   * @see MHD3DProjection
   */
  MHD3DProjectionPrimE(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~MHD3DProjectionPrimE();
  
  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Get extra variable names
   */
  virtual std::vector<std::string> getExtraVarNames() const;
  
  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal);
  
  /**
   * Set the matrix of the right and left eigenvectors
   * and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv,
					 RealVector& eValues,
					 const RealVector& normal);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  virtual void computePhysicalData(const Framework::State& state,
				   RealVector& data);
  
  /**
   * Set the total magnetic field and energy values
   */
  virtual void setDimensionalValuesPlusExtraValues(const Framework::State& state,
						   RealVector& result,
						   RealVector& extra);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data,
					    Framework::State& state);

protected:
  
  /// Computes the convective flux projected on a normal
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  /// Computes the physical convective flux
  virtual void computeStateFlux(const RealVector& pdata);
  
}; // end of class MHD3DProjectionPrimE
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Physics_MHD_MHD3DProjectionPrimE_hh
