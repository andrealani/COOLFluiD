#ifndef COOLFluiD_Physics_Poisson_PoissonConvVarSet_hh
#define COOLFluiD_Physics_Poisson_PoissonConvVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Poisson/PoissonConvTerm.hh"
#include "Framework/EquationSetData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for an Poisson physical model.
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 */
class PoissonConvVarSet : public Framework::ConvectiveVarSet {
public:
  typedef PoissonConvVarSet POISSONSET;
  typedef PoissonConvTerm PTERM;
  
  /**
   * Constructor
   * @see PoissonlModel
   */
  PoissonConvVarSet(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~PoissonConvVarSet();

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
   * Set the jacobians
   */
  virtual void computeJacobians()
  {
    throw Common::NotImplementedException (FromHere(),"PoissonConvVarSet::computeJacobians()");
  }
  
  /**
   * Compute the pressure derivative
   */
  virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp)
  {
    throw Common::NotImplementedException (FromHere(),"PoissonConvVarSet::computePressureDerivatives()");
  }
  
  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
                          const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"PoissonConvVarSet::splitJacobian()");
  }

  
  
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"PoissonConvVarSet::computeEigenValuesVectors()");
  }
  
  /// Set the PhysicalData corresponding to the given State
  virtual void computePhysicalData (const Framework::State& state, RealVector& pdata) = 0;
  
  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "PoissonConvVarSet::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "PoissonConvVarSet::setAdimensionalValues() not implemented");
  }

  /**
   * Get the normal speed
   */
  CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const
  {
    return 1.;
  }

  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);

  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);
  
  /**
   * Get the model
   */
  Common::SafePtr<PoissonConvTerm> getModel() const
  {
    return _model;
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
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
  
protected:  
  
  /// @returns the number of equations of this VarSet
  CFuint getNbEqs() const { return 1; std::cout << "getNbEqs()\n"; }
 
  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& pdata,
			   const RealVector& normals);

  /**
   * Compute the convective flux
   */
  virtual void computeStateFlux(const RealVector& vars);
  
private:
  
  /// acquaintance of the model
  Common::SafePtr<PoissonConvTerm> _model;
      
}; // end of class PoissonConvVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonConvVarSet_hh
