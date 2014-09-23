#ifndef COOLFluiD_Physics_NavierStokes_ArcJetLTEPhysicalModel_hh
#define COOLFluiD_Physics_NavierStokes_ArcJetLTEPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/NSTerm.hh"
#include "ArcJet/ArcJetReactionTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/ConvectionDiffusionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an ArcJetLTEPhysicalModel.
 *
 * @author Andrea Lani
 *
 */

template <int DIM>
class ArcJetLTEPhysicalModel :
	public Framework::ConvectionDiffusionReactionPM
	<Framework::MultiScalarTerm<NavierStokes::EulerTerm>,
	 NavierStokes::NSTerm, ArcJetReactionTerm<Framework::BaseTerm> > {
public:
  
   /**
    * Constructor without arguments
    */
  ArcJetLTEPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ArcJetLTEPhysicalModel();

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
    return std::string("ArcJetLTE" + Common::StringOps::to_str(DIM) + "D");
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

}; // end of class ArcJetLTEPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetLTEPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_ArcJetLTEPhysicalModel_hh
