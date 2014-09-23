#ifndef COOLFluiD_Physics_NEQ_Euler3DNEQCons_hh
#define COOLFluiD_Physics_NEQ_Euler3DNEQCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Physics {
    
    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 3D for conservative
   * variables with chemical NEQ
   *
   * @author Andrea Lani, Janos Molnar
   */
class Euler3DNEQCons : 
	public NavierStokes::MultiScalarVarSet<NavierStokes::Euler3DVarSet> {
public: // classes

  /**
   * Constructor
   * @see Euler3D
   */
  Euler3DNEQCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Euler3DNEQCons();

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
   * Set the jacobian matrix
   */
  void computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob);

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
  void computePhysicalData(const Framework::State& state,
			   RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data,
				    Framework::State& state);
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs);
  
private:

  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// constant of perfect gases
  CFreal _Rgas;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;

  /// array to store the mass fractions of elements
  RealVector _ys;
  
  /// temporary matrix of right eigenvalues
  RealMatrix _rightEv;
  
  /// temporary matrix of left eigenvalues
  RealMatrix _leftEv;
 
  /// useful coefficients
  RealVector _alpha;
  
  /// useful coefficients
  RealVector _RiGas;
  
  /// species molar masses
  RealVector _mmasses;
   
  /// f_i coefficients
  RealVector _fcoeff;

}; // end of class Euler3DNEQCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler3DNEQCons_hh
