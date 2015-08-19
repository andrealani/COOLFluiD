#ifndef COOLFluiD_Physics_Poisson_PoissonDiff3DCons_hh
#define COOLFluiD_Physics_Poisson_PoissonDiff3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Poisson/PoissonDiffVarSet.hh"
#include "Poisson/PoissonConvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Poisson physical model 3D for primitive
   * variables
   *
   * @author Alejandro Alvarez
   */
class PoissonDiff3DCons : public PoissonDiffVarSet {
public:

  typedef PoissonConvTerm PTERM;

  
  /**
   * Constructor
   * @see PoissonDiffVarSet
   */
  PoissonDiff3DCons(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~PoissonDiff3DCons();
  
  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  void setGradientVars(const std::vector<RealVector*>& states,
		       RealMatrix& values,
		       const CFuint stateSize);

  /**
   * Compute required gradients (velocity, Temperature) starting from the gradients of the states
   */
  void setGradientVarGradients(const std::vector<RealVector*>& states,
                               const std::vector< std::vector<RealVector*> >& stateGradients,
                               std::vector< std::vector<RealVector*> >& gradVarGradients,
                               const CFuint stateSize);

  /**
   * Compute the gradients of the states starting from gradient variable gradients (pressure, velocity, temperature)
   */
  void setStateGradients(const std::vector<RealVector*>& states,
                         const std::vector< std::vector<RealVector*> >& gradVarGradients,
                         std::vector< std::vector<RealVector*> >& stateGradients,
                         const CFuint stateSize);

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);
  
private:
  
  /// convective model
  Common::SafePtr<PTERM> m_eulerModel;

}; // end of class PoissonDiff3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonDiff3DCons_hh
