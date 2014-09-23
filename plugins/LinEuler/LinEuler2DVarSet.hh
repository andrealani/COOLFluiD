#ifndef COOLFluiD_Physics_LinEuler_LinEuler2DVarSet_hh
#define COOLFluiD_Physics_LinEuler_LinEuler2DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinEulerVarSet.hh"
#include "Framework/EquationSetData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 2D Lineariazed Euler physical model.
 * @author Lilla Edit Koloszar 
 * @author Nadege Villedieu
 * @author Tomas Kopacek
 * @author Matteo Parsani (modified for release 2009.3)
 */
class LinEuler2DVarSet : public LinEulerVarSet {
public: // classes

  typedef LinEuler2DVarSet LINEULERSET;

  /**
   * Constructor
   * @see LinEulerPhysicalModel
   */
  LinEuler2DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~LinEuler2DVarSet();

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
                          const RealVector& normal) = 0;

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal) = 0;

  /**
   * Set the first right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  virtual void setEigenVect1(RealVector& r1,
                             Framework::State& state,
                             const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"LinEuler2DVarSet::setEigenVect1()");
  }

  /**
   * Set the second right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  virtual void setEigenVect2(RealVector& r2,
                             Framework::State& state,
                             const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"LinEuler2DVarSet::setEigenVect2()");
  };

  /**
   * Set the third right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  virtual void setEigenVect3(RealVector& r3,
                             Framework::State& state,
                             const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"LinEuler2DVarSet::setEigenVect3()");
  }

  /**
   * Set the fourth right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  virtual void setEigenVect4(RealVector& r4,
                             Framework::State& state,
                             const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"LinEuler2DVarSet::setEigenVect4()");
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
  CFreal un = data[LinEulerTerm::u]*normal[0] + data[LinEulerTerm::v]*normal[1]; 
      return un;
  }

  /**
   * Get some data corresponding to the subset of equations related with
   * this variable set
   * @pre The most concrete ConvectiveVarSet will have to set these data
   */
  static std::vector<Framework::EquationSetData>& getEqSetData()
  {
    static std::vector<Framework::EquationSetData> eqSetData;
    return eqSetData;
  }

   /**
   * get the size of the extra physical vars
   */
  virtual CFuint getExtraPhysicalVarsSize()
  {
    return 4;
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(2); velIDs[XX] = 1; velIDs[YY] = 2; 
  }
  
protected:

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);

  /**
   * Compute the physical convective flux
   */
  virtual void computeStateFlux(const RealVector& pdata);

  /**
   * Get the number of equations of this VarSet
   */
  CFuint getNbEqs() const
  {
    return 4;
  }

  /**
   * Get the maximum eigen value
   */
  virtual CFreal getMaxEigenValue(const RealVector& pdata,
               const RealVector& normal);

  /**
   * Get the maximum absolute eigenvalue
   */
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata,
            const RealVector& normal);

  /**
   * Set the vector of the eigenValues
   */
  virtual void computeEigenValues(const RealVector& pdata,
               const RealVector& normal, RealVector& result);

}; // end of class LinEuler2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_LinEuler2DVarSet_hh
