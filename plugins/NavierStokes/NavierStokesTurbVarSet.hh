#ifndef COOLFluiD_Physics_NavierStokes_NavierStokesTurbVarSet_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokesTurbVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/MathConsts.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////
	
namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a diffusive variable set for the NavierStokes model
   * with turbulence
   *
   * @author Thomas Wuilbaut
   * @author Andrea Lani
   */
template <typename BASE, int SGROUP>
class NavierStokesTurbVarSet : public BASE {
public: // classes
  
  typedef NSTurbTerm DTERM; 
  
  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesTurbVarSet (const std::string& name, 
			  Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokesTurbVarSet();
  
  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Get the adimensional dynamic viscosity
   */
  virtual CFreal getLaminarDynViscosityFromGradientVars(const RealVector& state) = 0;
  
  /**
   * Get the adimensional turbulent dynamic viscosity
   */
  virtual CFreal getTurbDynViscosityFromGradientVars
  (const RealVector& state,
   const std::vector<RealVector*>& gradients) = 0;

  /**
   * Get the current adimensional dynamic viscosity
   */
  virtual CFreal getCurrDynViscosity()
  {
    return ((getModel().getPhysicalData())[NavierStokes::NSTurbTerm::MU] +
	    (getModel().getPhysicalData())[NavierStokes::NSTurbTerm::MUT]);
  }

  /// Compute the transport properties
  virtual void computeTransportProperties(const RealVector& state,
					  const std::vector<RealVector*>& gradients,
					  const RealVector& normal);
  
  /// Compute the stress tensor
  virtual void computeStressTensor(const RealVector& state,
				   const std::vector<RealVector*>& gradients,
				   const CFreal& radius); 
  
  /**
   * Get the model
   */
  NavierStokes::NSTurbTerm& getModel() {return *_turbModel;}
  
protected:
  
  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state) 
  {
    throw Common::NotImplementedException (FromHere(),"NavierStokesTurbVarSet::setGradientState()");
  }
  
protected:
  
  /// physical model
  Common::SafePtr<NavierStokes::NSTurbTerm> _turbModel;

  /// turbulent viscosity coefficient
  CFreal _turbDynViscCoeff;  
  
  /// variable ID (in the state) of theturbulent kinetic energy K
  CFuint _kID;
  
}; // end of class NavierStokesTurbVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes
} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesTurbVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokesTurbVarSet_hh
