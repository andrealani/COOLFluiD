#ifndef COOLFluiD_Physics_NavierStokes_Euler1DVarSet_hh
#define COOLFluiD_Physics_NavierStokes_Euler1DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "EulerVarSet.hh"
#include "Framework/EquationSetData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 1D Euler physical model.
 *
 * @author Andrea Lani
 */
class Euler1DVarSet : public EulerVarSet {
public: // classes

  typedef Euler1DVarSet EULERSET;

  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  Euler1DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Euler1DVarSet();

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
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const = 0;
  
  /**
   * Get the normal speed
   */
  CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const
  {
    return data[EulerTerm::VX]*normal[XX];
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
  
protected:
  
  /// Computes the convective flux projected on a normal
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  /// Computes the physical convective flux
  virtual void computeStateFlux(const RealVector& pdata);
  
  /// @returns the number of equations of this VarSet
  CFuint getNbEqs() const { return 3; }
  
}; // end of class Euler1DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler1DVarSet_hh
