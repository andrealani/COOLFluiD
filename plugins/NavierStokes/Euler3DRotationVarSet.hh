#ifndef COOLFluiD_Physics_NavierStokes_Euler3DRotationVarSet_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DRotationVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "EulerVarSet.hh"
#include "Framework/EquationSetData.hh"

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
class Euler3DRotationVarSet : public EulerVarSet {
public: // classes

  typedef Euler3DRotationVarSet EULERSET;

  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  Euler3DRotationVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Euler3DRotationVarSet();

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
   * Get the absolute speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DRotationVarSet::getSpeed()");
  }

  /**
   * Get the absolute normal speed
   */
  CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const
  {
    return data[EulerTerm::VX]*normal[XX] +
      data[EulerTerm::VY]*normal[YY] +
      data[EulerTerm::VZ]*normal[ZZ];
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

  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3;
  }
  
protected:

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& data,
			   const RealVector& normals);

  /**
   * Compute the convective flux
   */
  virtual void computeStateFlux(const RealVector& data);
  
  /**
   * Get the number of equations of this VarSet
   */
  CFuint getNbEqs() const
  {
    return 5;
  }

}; // end of class Euler3DRotationVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DRotationVarSet_hh
