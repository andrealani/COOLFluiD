#ifndef COOLFluiD_Physics_SA_NavierStokesSAPhysicalModel_hh
#define COOLFluiD_Physics_SA_NavierStokesSAPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/ConvectionDiffusionPM.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokesSAPhysicalModel.
 *
 * @author Andrea Lani
 *
 */

template <int DIM>
class NavierStokesSAPhysicalModel :
	public Framework::ConvectionDiffusionPM
	<Framework::MultiScalarTerm<NavierStokes::EulerTerm>,
	 NavierStokes::NSTurbTerm> {
public:

   /**
    * Constructor without arguments
    */
  NavierStokesSAPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesSAPhysicalModel();

   /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("NavierStokes" + Common::StringOps::to_str(DIM) + "DSA");
  }

  /**
   * Get the convective name
   */
  virtual  std::string getConvectiveName() const;

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const;

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

  /**
   * Set the reference value for time
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceTime();

}; // end of class NavierStokesSAPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

#include "NavierStokesSAPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_NavierStokesSAPhysicalModel_hh
