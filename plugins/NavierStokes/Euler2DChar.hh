#ifndef COOLFluiD_Physics_NavierStokes_Euler2DChar_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DChar_hh

//////////////////////////////////////////////////////////////////////////////

#include "Euler2DVarSet.hh"

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
class Euler2DChar : public Euler2DVarSet {
public:

  /**
   * Constructor
   * @see Euler3D
   */
  Euler2DChar(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Euler2DChar();

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;
  
  /// Set the PhysicalData corresponding to the given State
  virtual void computePhysicalData(const Framework::State& state, RealVector& data)
  {
    throw Common::NotImplementedException(FromHere(), "Euler2DChar::computePhysicalData() not implemented");
  }
  
  /// Set the State correspoding to the PhysicalData
  virtual void computeStateFromPhysicalData(const RealVector& pdata, Framework::State& state) 
  {
    throw Common::NotImplementedException(FromHere(), "Euler2DChar::computeStateFromPhysicalData() not implemented");
  }
  
  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Set the scalar part of the jacobian
   */
  void computeScalarJacobian(const RealVector& normal,
                      RealVector& jacob);

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
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const
  {
    throw Common::NotImplementedException
      (FromHere(), "Euler2DChar::getSpeed() not implemented");
  }

  /**
   * Give dimensional values to the adimensional state variables
   */
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "Euler2DChar::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  void setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "Euler2DChar::setAdimensionalValues() not implemented");
  }
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DChar::computePerturbedPhysicalData()");
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  { 
    velIDs.resize(2); velIDs[XX] = 1; velIDs[YY] = 2;
  }
  
}; // end of class Euler2DChar

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DChar_hh
