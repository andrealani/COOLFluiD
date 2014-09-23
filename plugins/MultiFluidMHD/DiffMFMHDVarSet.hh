#ifndef COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDVarSet_hh
#define COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarVarSetBase.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "MultiFluidMHD/DiffMFMHDTerm.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a diffusive variable set for the multi-fluid NavierStokes model
   *
   * @author Alejandro Alvarez
   * 
   */
class DiffMFMHDVarSet : public Framework::DiffusiveVarSet {
public: // classes

  typedef Framework::MultiScalarTerm<EulerMFMHDTerm> PTERM;
  
  /**
   * Constructor
   * @see NavierStokes
   */
  DiffMFMHDVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~DiffMFMHDVarSet();
  
  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  virtual void setComposition(const RealVector& state,
                              const bool isPerturb,
                              const CFuint iVar) = 0;
  
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
   * Get the dimensional dynamic viscosity
   */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients) = 0;
  
  /**
   * Get the current adimensional dynamic viscosity
   */
  virtual CFreal getCurrDynViscosity()
  {
    const CFuint nbSpecies = getModel().getNbSpecies();
    CFreal currDynViscosity = 0.;
    for (CFuint i = 0; i < nbSpecies; ++i) {
      currDynViscosity = std::max(currDynViscosity, (getModel().getPhysicalData())[i]);
    }  
    return currDynViscosity;
  }
  
  /**
   * Get the current adimensional thermal conductivity
   */
  virtual CFreal getCurrThermConductivity()
  {
    throw Common::NotImplementedException(FromHere(),"getCurrThermConductivity");
    //     return (getModel().getPhysicalData())[DiffMFMHDTerm::LAMBDA];
  }
  
  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS)
  {
    return (iEqSS == 1) ? getCurrDynViscosity() : 0.;
  }
  
  /**
   * Get the heat flux
   */
  virtual RealVector& getHeatFlux(const RealVector& state,
                             const std::vector<RealVector*>& gradients,
                             const RealVector& normal);
  
  /**
   * Get the adimensional density
   */
  virtual CFreal getDensity(const RealVector& state) = 0;
  
  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
                                      const CFreal& dynViscosity) = 0;
  
  /**
   * Get the model
   */
  DiffMFMHDTerm& getModel() {return *m_model;}
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius) = 0;
  
  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius) = 0;
  
  /// Compute the stress tensor
  virtual void computeStressTensor(const RealVector& state,
				   const std::vector<RealVector*>& gradients,
				   const CFreal& radius); 
  /**
   * Set the wall distance
   */
  void setWallDistance(CFreal dist) {_wallDistance = dist;}
  
protected:
  
  /// Compute the transport properties
  virtual void computeTransportProperties(const RealVector& state,
					  const std::vector<RealVector*>& gradients,
					  const RealVector& normal);
  
  /// Compute the braginskii thermal conductivity coefficients
  void computeBraginskiiThermConduct(const RealVector& state);
  
  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state) = 0;
 
  /**
   * Get the dimensional dynamic viscosity
   */
  virtual RealVector& getDynViscosityVec(const RealVector& state, const std::vector<RealVector*>& gradients) = 0;
  
 
  /**
   * Get the dimensional dynamic viscosity
   */
  virtual RealVector& getThermConductivityVec(const RealVector& state,
			      const CFreal& dynViscosity) = 0;
  
private:

  /// physical model
  Common::SafePtr<DiffMFMHDTerm> m_model;
  
  /// convective model
  Common::SafePtr<PTERM> m_eulerModel;

protected:
  
  /// 2./3.
  const CFreal _twoThird;
  
  /// flag to use back up values for certain quantities
  bool _useBackUpValues;

  /// flag to set back up values for certain quantities
  bool _setBackUpValues;
  
  /// wall distance
  CFreal _wallDistance;
  
  /// gradient variables
  RealVector _gradState;
  
  /// local normal
  RealVector _normal;
  
  /// Heat Flux vector of the species
  RealVector _heatFlux;
  
  /// back up value for dynamic viscosity coefficient
  RealVector m_dynViscCoeff;

  /// back up value for thermal conductivity coefficient
  RealVector m_thermCondCoeff;
    
  /// gradient variables
  std::vector<RealMatrix> m_tau;
  
  ///Curr Dynamic Viscosity
  RealVector m_currDynViscosity;
  
  /// variable ID (in the state) of the Vx component of the velocity
  RealVector m_uID;
  
  /// variable ID (in the state) of the Vy component of the velocity
  RealVector m_vID;

  /// variable ID (in the state) of the Vz component of the velocity
  RealVector m_wID;
  
  /// variable ID (in the state) of the temperature
  RealVector m_TID;
  
  /// lambda coefficient for heat Flux computation
  RealVector m_lambda;
  
  ///Thermal conductivity parallel to the magnetic field (Braginskii)
  CFreal _kappaParallel;
  
  ///Thermal conductivity perpendicular to the magnetic field (Braginskii) 
  CFreal _kappaPerpendicular;
  
  ///Thermal conductivity perpendicular to the magnetic field (Braginskii) 
  CFreal _betaWedge;  
  
  
}; // end of class DiffMFMHDVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDVarSet_hh
