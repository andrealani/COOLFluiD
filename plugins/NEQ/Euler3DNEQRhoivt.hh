#ifndef COOLFluiD_Physics_NEQ_Euler3DNEQRhoivt_hh
#define COOLFluiD_Physics_NEQ_Euler3DNEQRhoivt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 2D for primitive
   * variables with chemical NEQ
   *
   * @author Andrea Lani
   */
class Euler3DNEQRhoivt : 
	public NavierStokes::MultiScalarVarSet<NavierStokes::Euler3DVarSet> {
public: // classes

  /**
   * Constructor
   * @see Euler3D
   */
  Euler3DNEQRhoivt(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Euler3DNEQRhoivt();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Get extra variable names
   */
  virtual std::vector<std::string> getExtraVarNames() const;

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;
  
  /**
   * Set the jacobians
   */
  void computeJacobians();
 
  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			  RealMatrix& jacobMin,
			  RealVector& eValues,
			  const RealVector& normal);
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal);
  
  /**
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const;

  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
				    RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
				     RealVector& result);

  /**
   * Set other adimensional values for useful physical quantities
   */
  virtual void setDimensionalValuesPlusExtraValues
  (const Framework::State& state, RealVector& result,
   RealVector& extra);
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  virtual void computePhysicalData(const Framework::State& state,
				   RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data,
					    Framework::State& state);
  /**
   * Compute the pressure derivative
   */
  virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp);
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs);
  
protected:
  
  /**
   * Set all thermodynamic quantities in the physical data array
   */
  virtual void setThermodynamics(CFreal rho,
				 const Framework::State& state,
				 RealVector& data);
  
  /**
   * Get the ID of the temperature given the number of species
   */
  CFuint getTempID(CFuint nbSpecies) const
  {
    return nbSpecies + 3;
  }
  
protected:
  
  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;
  
  /// array to store the mass fractions of elements
  RealVector _ye;
  
  /// array to store Rgas/molar mass for each species
  RealVector _Rspecies;

}; // end of class Euler3DNEQRhoivt

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler3DNEQRhoivt_hh
