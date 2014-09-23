#ifndef COOLFluiD_Physics_NavierStokes_Euler3DCons_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 3D for conservative
 * variables
 *
 * @author Andrea Lani
 */
class Euler3DCons : public Euler3DVarSet {
public: // classes

public: // methods

  /**
   * Constructor
   * @see Euler3D
   */
  Euler3DCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Euler3DCons();

  /**
   * Set up private data and data
   */
  virtual void setup();

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();
  
  /**
   * Set the jacobian matrix
   */
  virtual void computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob);
  
  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal);
  
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv,
					 RealVector& eValues,
					 const RealVector& normal);
  
  /**
   * Set the first right eigen vector (corresponding to \f$\vec{u}\cdot\vec{n}\f$}
   */
  void setEigenVect1(RealVector& r1,
		     Framework::State& state,
                     const RealVector& normal);
  /**
   * Set the second right eige�n vector (corresponding to \f$\vec{u}\cdot\vec{n}\f$}
   */
  void setEigenVect2(RealVector& r2,
                     Framework::State& state,
                     const RealVector& normal);
  /**
   * Set the third right eigen vector (corresponding to \f$\vec{u}\cdot\vec{n}\f$}
   */
  void setEigenVect3(RealVector& r3,
                     Framework::State& state,
                     const RealVector& normal);

  /**
   * Set the fourth right eigen vector (corresponding to \f$\vec{u}\cdot\vec{n}+ a \f$}
   */
  void setEigenVect4(RealVector& r4,
                     Framework::State& state,
                     const RealVector& normal);

  /**
   * Set the fifth right eigen vector (corresponding to \f$\vec{u}\cdot\vec{n}- a \f$}
   */
  void setEigenVect5(RealVector& r5,
                     Framework::State& state,
                     const RealVector& normal);

  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
				    RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                             RealVector& result);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  virtual void computePhysicalData(const Framework::State& state,
			   RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data,
				    Framework::State& state);

  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DCons::computePerturbedPhysicalData()");
  }
  
  /**
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const;
   
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3;
  }
  
private:

  /**
   * Set the constant part (independent form the solution) of the
   * jacobians
   */
  virtual void setConstJacob();
  
private:

  /// temporary matrix of right eigenvalues
  RealMatrix                       _rightEv;

  /// temporary matrix of left eigenvalues
  RealMatrix                       _leftEv;

}; // end of class Euler3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DCons_hh
