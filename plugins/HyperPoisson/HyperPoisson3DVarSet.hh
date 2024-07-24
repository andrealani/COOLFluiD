#ifndef COOLFluiD_Physics_HyperPoisson_HyperPoisson3DVarSet_hh
#define COOLFluiD_Physics_HyperPoisson_HyperPoisson3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "HyperPoissonVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 3D Hyperbolized Poisson physical model.
 *
 * @author Rayan Dhib
 * @author Andrea Lani
 */
class HyperPoisson3DVarSet : public HyperPoissonVarSet {
public: // classes
  
  /**
   * Constructor
   * @see HyperPoissonPhysicalModel
   */
  HyperPoisson3DVarSet(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~HyperPoisson3DVarSet();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const = 0;

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal) {} 
  
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv,
					 RealVector& eValues,
					 const RealVector& normal) {} 
  
  /**
   * Set the first right eigen vector
   */
  virtual void setEigenVect1(RealVector& r1,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"HyperPoisson3DVarSet::setEigenVect1()");
  }
  
  /**
   * Set the second right eigen vector
   */
  virtual void setEigenVect2(RealVector& r2,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"HyperPoisson3DVarSet::setEigenVect2()");
  };

  /**
   * Set the third right eigen vector
   */
  virtual void setEigenVect3(RealVector& r3,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"HyperPoisson3DVarSet::setEigenVect3()");
  }
  
  /**
   * Set the fourth right eigen vector
   */
  virtual void setEigenVect4(RealVector& r4,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"HyperPoisson3DVarSet::setEigenVect4()");
  }
  
  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);
  
protected:

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& pdata,
			   const RealVector& normals);

  /**
   * Compute the convective flux
   */
  virtual void computeStateFlux(const RealVector& vars);
  
  /**
   * Get the number of equations of this VarSet
   */
  CFuint getNbEqs() const {return 4;}
  
}; // end of class HyperPoisson3DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoissonVarSet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_HyperPoisson_HyperPoisson3DVarSet_hh
