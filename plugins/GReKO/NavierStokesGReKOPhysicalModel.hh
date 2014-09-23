#ifndef COOLFluiD_Physics_GReKO_NavierStokesGReKOPhysicalModel_hh
#define COOLFluiD_Physics_GReKO_NavierStokesGReKOPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NSTurbTerm.hh"
#include "GReKOReactionTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/ConvectionDiffusionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokesGReKOPhysicalModel.
 *
 * @author Khalil Bensassi
 *
 */

template <int DIM>
class NavierStokesGReKOPhysicalModel :
	public Framework::ConvectionDiffusionReactionPM
	<Framework::MultiScalarTerm<NavierStokes::EulerTerm>,
	 NavierStokes::NSTurbTerm, GReKOReactionTerm> {
public:

   /**
    * Constructor without arguments
    */
  NavierStokesGReKOPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesGReKOPhysicalModel();

   /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure (Config::ConfigArgs& args);

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("NavierStokes" +  Common::StringOps::to_str(DIM) + "DGReKO");
  }

  /**
   * Get the convective name
   */
  virtual  std::string getConvectiveName() const;

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const
  {
    return getTypeName();
  }

  /**
   * Get the source name
   */
  virtual std::string getSourceName() const
  {
    return getTypeName();
  }

  /**
   * @return the space dimension of the SubSystem
   */
  virtual CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  virtual CFuint getNbEquations() const;

  /**
   * Set the reference values
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceValues();

}; // end of class NavierStokesGReKOPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace GReKO

  } // end of namespace Physics

} // end of namespace COOLFluiD

#include "NavierStokesGReKOPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokesGReKOPhysicalModel_hh
