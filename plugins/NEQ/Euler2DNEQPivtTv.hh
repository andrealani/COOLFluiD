#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQPivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQPivtTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "NEQ/Euler2DNEQPivt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 2D for primitive
   * variables with chemical and thermal NEQ
   *
   * @author Andrea Lani
   */
class Euler2DNEQPivtTv : public Euler2DNEQPivt {
  
public: // classes
  
  /**
   * Constructor
   * @see EulerNEQ
   */
  Euler2DNEQPivtTv(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~Euler2DNEQPivtTv();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();
  
  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;
  
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
  virtual void computeEigenValuesVectors(RealMatrix& rightTv,
				     RealMatrix& leftTv,
				     RealVector& eValues,
				     const RealVector& normal);
  
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
  /**
   * Checks validity of data
   * @pre data is assumed to be of number of equations size
   */
   bool isValid(const RealVector& data);
    
protected:
  
  /**
   * Set enthalpy, energy, sound speed taking into account
   * thermal non-equilibrium effects
   */
  virtual void setThermodynamics(CFreal rho, 
				 const Framework::State& state, 
				 RealVector& data);
  
  /// get the electron temperature
  virtual CFreal getTe(const Framework::State& state); 
  
protected:
  
  /// array with all different vibrational dimensional temperatures
  RealVector _tvDim;
  
  /// molecules iDs
  std::vector<CFuint> _moleculesIDs;
  
}; // end of class Euler2DNEQPivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQPivtTv_hh
