#ifndef COOLFluiD_Physics_LinEuler_LinEuler3DCons_hh
#define COOLFluiD_Physics_LinEuler_LinEuler3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinEuler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;
  }

  namespace Physics {

    namespace  LinearizedEuler{


//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Linearized Euler physical model 3D for conservative
 * variables
 * @author Kris Van den Abeele
 * @author Matteo Parsani (modified for release 2009.3)
 */
class LinEuler3DCons : public LinEuler3DVarSet {
public: // classes

  /**
   * Constructor
   * @see LinEuler2D
   */
  LinEuler3DCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~LinEuler3DCons();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
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
   * Set the jacobian matrix
   */
  void computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob);

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
   * Set the first right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  void setEigenVect1(RealVector& r1, Framework::State& state, const RealVector& normal);
  /**
   * Set the second right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  void setEigenVect2(RealVector& r2, Framework::State& state, const RealVector& normal);

  /**
   * Set the third right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n} \f$)
   */
  void setEigenVect3(RealVector& r3, Framework::State& state, const RealVector& normal);

  /**
   * Set the fourth right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}- a \f$)
   */
  void setEigenVect4(RealVector& r4, Framework::State& state, const RealVector& normal);

  /**
   * Set the fifth right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  void setEigenVect5(RealVector& r5, Framework::State& state, const RealVector& normal);

  /**
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const;

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
   * Compute the perturbed states data
   */
  void computePerturbedStatesData(const std::vector<Framework::State*>& states,
                                  const CFuint nbStatesInVec,
                                  const CFuint iVar);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see LinEulerPhysicalModel
   */
  void computePhysicalData(const Framework::State& state, RealVector& data);

  /**
   * Set a State starting from the given PhysicalData
   * @see LinEulerPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data, Framework::State& state);

protected:

  /**
   * Set the constant part (independent form the solution) of the
   * jacobians
   */
  void setConstJacob();
  
//   virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
//   
//   virtual void computeStateFlux(const RealVector& pdata);

private: // data

  /// temporary matrix of right eigenvalues
  RealMatrix                       _rightEv;

  /// temporary matrix of left eigenvalues
  RealMatrix                       _leftEv;

}; // end of class LinEuler3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_LinEuler3DCons_hh
