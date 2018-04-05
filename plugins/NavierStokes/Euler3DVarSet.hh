#ifndef COOLFluiD_Physics_NavierStokes_Euler3DVarSet_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 3D Euler physical model.
 *
 * @author Andrea Lani
 */
class Euler3DVarSet : public EulerVarSet {
public: // classes
  
  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  Euler3DVarSet(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~Euler3DVarSet();

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
   * Set the first right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  virtual void setEigenVect1(RealVector& r1,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DVarSet::setEigenVect1()");
  }
  
  /**
   * Set the second right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  virtual void setEigenVect2(RealVector& r2,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DVarSet::setEigenVect2()");
  };

  /**
   * Set the third right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  virtual void setEigenVect3(RealVector& r3,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DVarSet::setEigenVect3()");
  }
  
  /**
   * Set the fourth right eigen vector (corresponding to \f$\vec{u}\cdot\vec{n}+ a \f$}
   */
  virtual void setEigenVect4(RealVector& r4,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DVarSet::setEigenVect4()");
  }

  /**
   * Set the fifth right eigen vector (corresponding to \f$\vec{u}\cdot\vec{n}- a \f$}
   */
  virtual void setEigenVect5(RealVector& r5,
			     Framework::State& state,
			     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DVarSet::setEigenVect5()");
  }

  /**
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DVarSet::getSpeed()");
  }

  /**
   * Get the normal speed
   */
  CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const
  {
    if (normal.size() == 3) {
      return data[EulerTerm::VX]*normal[XX] +
	data[EulerTerm::VY]*normal[YY] +
	data[EulerTerm::VZ]*normal[ZZ];
    }
    // 2D and 1/2 
    return data[EulerTerm::VX]*normal[XX] + data[EulerTerm::VY]*normal[YY];
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
  CFuint getNbEqs() const {return 5;}
  
}; // end of class Euler3DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DVarSet_hh
