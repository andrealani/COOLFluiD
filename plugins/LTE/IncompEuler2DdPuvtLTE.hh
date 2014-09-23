#ifndef COOLFluiD_Physics_LTE_IncompEuler2DdPuvtLTE_hh
#define COOLFluiD_Physics_LTE_IncompEuler2DdPuvtLTE_hh

//////////////////////////////////////////////////////////////////////////////

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
   * This class represents an incompresible Euler physical model 
   * 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class IncompEuler2DdPuvtLTE : public NavierStokes::IncompEuler2DVarSet {
public: // classes

  /**
   * Constructor
   * @see IncompEuler2D
   */
  IncompEuler2DdPuvtLTE(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~IncompEuler2DdPuvtLTE();

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
   * Set other adimensional values for useful physical quantities
   */
  virtual void setDimensionalValuesPlusExtraValues
  (const Framework::State& state, RealVector& result,
   RealVector& extra);

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;
  
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
   * Compute the perturbed states data
   */
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
  
private:

  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// array to store density, enthalpy and energy
  RealVector _dhe;

  /// array to store the volume composition for each species
  RealVector _x;

}; // end of class IncompEuler2DdPuvtLTE

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_IncompEuler2DdPuvtLTE_hh
