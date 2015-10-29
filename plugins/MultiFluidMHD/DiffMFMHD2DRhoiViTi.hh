#ifndef COOLFluiD_Physics_MultiFluidMHD_DiffMFMHD2DRhoiViTi_hh
#define COOLFluiD_Physics_MultiFluidMHD_DiffMFMHD2DRhoiViTi_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarVarSetBase.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a MultiFluidMHD Diffusive physical model 3D for primitive
   * variables
   *
   *
   * @author Alejandro Alvarez
   */
class DiffMFMHD2DRhoiViTi : public DiffMFMHD2DVarSet {
public:

  typedef Framework::MultiScalarTerm<EulerMFMHDTerm> PTERM;

  
  /**
   * Constructor
   * @see DiffMFMHD2D
   */
  DiffMFMHD2DRhoiViTi(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~DiffMFMHD2DRhoiViTi();

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  void setComposition(const RealVector& state,
		      const bool isPerturb,
		      const CFuint iVar)
  {
  }
  
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

  /**
   * Get the adimensional dynamic viscosity
   */
  CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients);
  
    /**
   * Get the adimensional thermal conductivity
   */
  CFreal getThermConductivity(const RealVector& state,const CFreal& dynViscosity);
  
  /**
   * Get the adimensional density
   */
  CFreal getDensity(const RealVector& state);

  /**
   * Get the adimensional thermal conductivity
   */
  RealVector& getThermConductivityVec(const RealVector& state,
			      const CFreal& dynViscosity)
  {
          return getModel().getThermConductivityDim();
  }
  
    /**
   * Get the adimensional dynamic viscosity
   */
  RealVector& getDynViscosityVec(const RealVector& state, const std::vector<RealVector*>& gradients);

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);
  
private:
  
  /// convective model
  Common::SafePtr<PTERM> m_eulerModel;
  
  RealVector m_density;
  
  RealVector m_m_i;

}; // end of class DiffMFMHD2DRhoiViTi

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_DiffMFMHD2DRhoiViTi_hh
