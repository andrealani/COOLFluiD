#ifndef COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHDVarSet_hh
#define COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHDVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarVarSetBase.hh"
#include "EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for an MultiFluidMHD physical model.
 *
 * @author Andrea Lani
 */
template <typename BASE>      
class MultiFluidMHDVarSet : public BASE {
public: // classes
  
  typedef Framework::MultiScalarTerm<EulerMFMHDTerm> PTERM;  
  
  /**
   * Constructor
   * @see MultiFluidMHDModel
   */
  MultiFluidMHDVarSet(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~MultiFluidMHDVarSet();

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
    throw Common::NotImplementedException (FromHere(),"MultiFluidMHDVarSet::computeJacobians()");
  }
  
  /**
   * Compute the pressure derivative
   */
  virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp)
  {
    throw Common::NotImplementedException (FromHere(),"MultiFluidMHDVarSet::computePressureDerivatives()");
  }
  
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
					 const RealVector& normal);
  
  /// Set the PhysicalData corresponding to the given State
  virtual void computePhysicalData (const Framework::State& state, RealVector& pdata) = 0;
  
  /// Get the speed
  virtual CFreal getSpeed(const Framework::State& state) const = 0;
  
  /// Get the normal speed
  CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const;
  
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
  
  /// Give dimensional values to the adimensional state variables
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result) = 0;
  /// Give adimensional values to the dimensional state variables
  virtual void setAdimensionalValues(const Framework::State& state,
				     RealVector& result) = 0;
				     
  
  /// Get the model
  Common::SafePtr<PTERM> getModel() const
  {
    cf_assert(m_mfModel.isNotNull());
    return m_mfModel;
  }
  
protected:
  
  /// Computes the convective flux projected on a normal
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  /// Computes the physical convective flux
  virtual void computeStateFlux(const RealVector& pdata);
  
private:
  
  /// acquaintance of the model
  Common::SafePtr<PTERM> m_mfModel;
      
}; // end of class MultiFluidMHDVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MultiFluidMHDVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHDVarSet_hh
