#ifndef COOLFluiD_Physics_NavierStokes_NavierStokesVarSet_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokesVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "NavierStokes/NSTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a diffusive variable set for the NavierStokes model
   *
   * @author Andrea Lani
   */
class NavierStokesVarSet : public Framework::DiffusiveVarSet {
public: // classes
  
  typedef NSTerm DTERM;
  
  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesVarSet();
  
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
   * Get the adimensional dynamic viscosity
   */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients) = 0;
  
  /**
   * Get the current adimensional dynamic viscosity
   */
  virtual CFreal getCurrDynViscosity()
  {
    return (getModel().getPhysicalData())[NSTerm::MU];
  }
  
  /**
   * Get the current adimensional thermal conductivity
   */
  virtual CFreal getCurrThermConductivity()
  {
    return (getModel().getPhysicalData())[NSTerm::LAMBDA];
  }
  
  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS)
  {
    return (getModel().getPhysicalData())[NSTerm::MU];
  }
  
  /**
   * Get the heat flux
   */
  virtual CFreal getHeatFlux(const RealVector& state,
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
  NSTerm& getModel() {return *_model;}
  
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
  
  /**
   * Get the axisymmetric source term
   */
  virtual void getAxiSourceTerm(const RealVector& physicalData,
				const RealVector& state,
				const std::vector<RealVector*>& gradients,
				const CFreal& radius,
				RealVector& source)
  {
    throw Common::NotImplementedException (FromHere(),"NavierStokesVarSet::getAxiSourceTerm()");
  }

protected:
  
  /// Compute the transport properties
  virtual void computeTransportProperties(const RealVector& state,
					  const std::vector<RealVector*>& gradients,
					  const RealVector& normal);
  
  /**
   * Set the flag asking to use back up values for certain quantities
   * instead of recomputing them
   */
  void useBackUpValues(bool flag) {_useBackUpValues = flag;}

  /**
   * Sset the flag asking to set back up values for certain quantities
   */
  void setBackUpValues(bool flag) {_setBackUpValues = flag;}
  
  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state) = 0;

private:

  /// physical model
  Common::SafePtr<NSTerm> _model;

protected:

  /// 2./3.
  const CFreal _twoThird;

  /// back up value for dynamic viscosity coefficient
  CFreal _dynViscCoeff;

  /// back up value for thermal conductivity coefficient
  CFreal _thermCondCoeff;

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
  
  /// gradient variables
  RealMatrix _tau;

  /// variable ID (in the state) of the Vx component of the velocity
  CFuint _uID;

  /// variable ID (in the state) of the Vy component of the velocity
  CFuint _vID;

  /// variable ID (in the state) of the Vz component of the velocity
  CFuint _wID;

  /// variable ID (in the state) of the temperature
  CFuint _TID;
  
}; // end of class NavierStokesVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokesVarSet_hh
