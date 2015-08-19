#ifndef COOLFluiD_Physics_Poisson_PoissonDiffVarSet_hh
#define COOLFluiD_Physics_Poisson_PoissonDiffVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Poisson/PoissonDiffTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework  {
   class PhysicalChemicalLibrary;
  }

  namespace Physics {
    
    namespace Poisson {
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Poisson physical model 2D for primitive
   * variables
   *
   * @author Alejandro Alvarez
   */

class PoissonDiffVarSet : public Framework::DiffusiveVarSet {
public: // classes

  /**
   * Constructor
   */
  PoissonDiffVarSet(const std::string& name,
		   Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */

  virtual ~PoissonDiffVarSet();
  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();
  
  
  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  virtual void setGradientVars(const std::vector<RealVector*>& states,
                               RealMatrix& values,
                               const CFuint stateSize) = 0;

  /**
   * Compute required gradients (pressure, velocity, temperature) starting from the gradients of the states
   */
  virtual void setGradientVarGradients(const std::vector<RealVector*>& states,
                                       const std::vector< std::vector<RealVector*> >& stateGradients,
                                       std::vector< std::vector<RealVector*> >& gradVarGradients,
                                       const CFuint stateSize) = 0;

  /**
   * Compute the gradients of the states starting from gradient variable gradients (pressure, velocity, temperature)
   */
  virtual void setStateGradients(const std::vector<RealVector*>& states,
                                 const std::vector< std::vector<RealVector*> >& gradVarGradients,
                                 std::vector< std::vector<RealVector*> >& stateGradients,
                                 const CFuint stateSize) = 0;
  

  /**
   * Get the model
   */
  PoissonDiffTerm& getModel() {return *_model;}
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal);


  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients)
  {
    throw Common::NotImplementedException (FromHere(),"PoissonDiffVarSet::getFlux()");
  }


  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);

protected:
  /// gradient variables
  RealVector _gradState;

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state) = 0;


private:

  /// pointer to the physical chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// physical model
  Common::SafePtr<PoissonDiffTerm> _model;
   
}; // end of class PoissonDiffVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonDiffVarSet_hh
