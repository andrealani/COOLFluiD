#ifndef COOLFluiD_Physics_IncompLTE_Euler2DLTEDemixCons_hh
#define COOLFluiD_Physics_LTE_IncompEuler2DLTEDemixCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/IncompEuler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Physics {
    
    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 2D for primitive
   * variables with LTE Demixing
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
class IncompEuler2DLTEDemixCons : public 
NavierStokes::MultiScalarVarSet<NavierStokes::IncompEuler2DVarSet> {
public: // classes

  /**
   * Constructor
   * @see Euler2D
   */
  IncompEuler2DLTEDemixCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~IncompEuler2DLTEDemixCons();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  void setup();

  /**
   * Get extra variable names
   */
  std::vector<std::string> getExtraVarNames() const;

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

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
  CFreal getSpeed(const Framework::State& state) const;

  /**
   * Give dimensional values to the adimensional state variables
   */
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  void setAdimensionalValues(const Framework::State& state,
                             RealVector& result);

  /**
   * Set other adimensional values for useful physical quantities
   */
  void setDimensionalValuesPlusExtraValues
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
    
private:

 /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// array to store density, enthalpy and energy
  RealVector _dhe;

  /// array to store the volume composition for each species
  RealVector _x;

  /// array to store the mass fractions of elements
  RealVector _ye;

}; // end of class Euler2DLTEDemixCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_IncompEuler2DLTEDemixCons_hh
