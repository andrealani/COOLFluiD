#ifndef COOLFluiD_Physics_NavierStokes_ICPLTEPhysicalModel_hh
#define COOLFluiD_Physics_NavierStokes_ICPLTEPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/NSTerm.hh"
#include "ICP/ICPReactionTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/ConvectionDiffusionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an ICPLTEPhysicalModel.
 *
 * @author Radek Honzatko
 *
 */

template <int DIM>
class ICPLTEPhysicalModel :
	public Framework::ConvectionDiffusionReactionPM
	<Framework::MultiScalarTerm<NavierStokes::EulerTerm>,
	 NavierStokes::NSTerm, ICPReactionTerm<Framework::BaseTerm> > {
public:
  
   /**
    * Constructor without arguments
    */
  ICPLTEPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ICPLTEPhysicalModel();

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  virtual std::string getTypeName() const
  {
    return std::string("ICPLTE" + Common::StringOps::to_str(DIM) + "D");
  }

  /**
   * Get the convective name
   */
  std::string getConvectiveName() const;

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

 private:

  /**
   * Set the reference values
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceValues();

}; // end of class ICPLTEPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICPLTEPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_ICPLTEPhysicalModel_hh
