#ifndef COOLFluiD_Physics_NavierStokes_Euler2DVarSet_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 2D Euler physical model.
 *
 * @author Andrea Lani
 */
class Euler2DVarSet : public EulerVarSet {
public: // classes

  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  Euler2DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Euler2DVarSet();

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
    throw Common::NotImplementedException (FromHere(),"Euler2DVarSet::setEigenVect1()");
  }

  /**
   * Set the second right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  virtual void setEigenVect2(RealVector& r2,
                             Framework::State& state,
                             const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DVarSet::setEigenVect2()");
  };

  /**
   * Set the third right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  virtual void setEigenVect3(RealVector& r3,
                             Framework::State& state,
                             const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DVarSet::setEigenVect3()");
  }

  /**
   * Set the fourth right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  virtual void setEigenVect4(RealVector& r4,
                             Framework::State& state,
                             const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DVarSet::setEigenVect4()");
  }
  
  /**
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const = 0;

  /**
   * Get the normal speed
   */
  CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const
  {
    return data[EulerTerm::VX]*normal[XX] + data[EulerTerm::VY]*normal[YY];
  }

  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);
  
protected:

  /// Computes the convective flux projected on a normal
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  /// Computes the physical convective flux
  virtual void computeStateFlux(const RealVector& pdata);
  
  /// @returns the number of equations of this VarSet
  CFuint getNbEqs() const { return 4; }
  
}; // end of class Euler2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DVarSet_hh
