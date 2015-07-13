#ifndef COOLFluiD_Physics_SA_NavierStokesSAPvtLTE_hh
#define COOLFluiD_Physics_SA_NavierStokesSAPvtLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "SA/NavierStokesSAPvt.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model for primitive
   * variables and SA turbulence model with LTE assumption
   *
   * @author Andrea Lani
   */
template <typename BASE>      
class NavierStokesSAPvtLTE : public NavierStokesSAPvt<BASE> {
public: // classes
  
  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesSAPvtLTE(const std::string& name,
		       Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesSAPvtLTE();
  
  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  virtual void setComposition(const RealVector& state,
			      const bool isPerturb,
			      const CFuint iVar);
  
  /**
   * Get the adimensional laminar dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getLaminarDynViscosityFromGradientVars(const RealVector& state);
  
  /**
   * Get the adimensional density
   * @pre the composition will be set here
   */
  virtual CFreal getDensity(const RealVector& state);
  
  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
				      const CFreal& dynViscosity);
  
  /**
   * set up private data
   */
  virtual void setup();
  
protected:
  
  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// array for the composition
  RealVector m_tempX;
  
}; // end of class NavierStokesSAPvtLTE
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesSAPvtLTE.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_NavierStokesSAPvtLTE_hh
