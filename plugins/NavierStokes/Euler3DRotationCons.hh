#ifndef COOLFluiD_Physics_NavierStokes_Euler3DRotationCons_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DRotationCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/Euler3DRotationVarSet.hh"

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
class Euler3DRotationCons : public Euler3DRotationVarSet {
public: // classes

public: // methods

  /**
   * Constructor
   * @see Euler3D
   */
  Euler3DRotationCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Euler3DRotationCons();

  /**
   * Set up private data and data
   */
  virtual void setup();

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
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  void setAdimensionalValues(const Framework::State& state,
                             RealVector& result);

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

  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DRotationCons::computePerturbedPhysicalData()");
  }
  
private:

  /**
   * Set the constant part (independent form the solution) of the
   * jacobians
   */
  void setConstJacob();
  
private:

  /// temporary matrix of right eigenvalues
  RealMatrix                       _rightEv;

  /// temporary matrix of left eigenvalues
  RealMatrix                       _leftEv;

}; // end of class Euler3DRotationCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DRotationCons_hh
