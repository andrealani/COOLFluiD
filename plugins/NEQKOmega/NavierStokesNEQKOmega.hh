#ifndef COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmega_hh
#define COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmega_hh

//////////////////////////////////////////////////////////////////////////////

#include "NEQ/NavierStokesNEQ.hh"
#include "KOmega/KOmegaReactionTerm.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokes NEQ KOmega model.
 *
 * @author Andrea Lani
 *
 */
template <int DIM, typename RTERM>
class NavierStokesNEQKOmega : 
	public NEQ::NavierStokesNEQ<DIM, 
				    Framework::MultiScalarTerm<NavierStokes::EulerTerm>, 
				    NavierStokes::NSTurbTerm,
				    KOmega::KOmegaReactionTerm<RTERM> > {
  
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  NavierStokesNEQKOmega(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesNEQKOmega();
  
  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  virtual std::string getTypeName() const
  {
    return std::string("NavierStokes" +  Common::StringOps::to_str(DIM) + "D" + "NEQKOmega");
  }

  /**
   * Get the convective name
   */
  virtual std::string getConvectiveName() const;
  
  /**
   * @return the number of equations of the SubSystem
   */
  virtual CFuint getNbEquations() const;
  
protected:

  /**
   * Set the reference values
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceValues();

}; // end of class NavierStokesNEQKOmega

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQKOmega.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmega_hh
