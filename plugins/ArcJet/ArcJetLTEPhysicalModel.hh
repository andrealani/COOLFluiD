#ifndef COOLFluiD_Physics_NavierStokes_ArcJetLTEPhysicalModel_hh
#define COOLFluiD_Physics_NavierStokes_ArcJetLTEPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "ArcJet/ArcJetTerm.hh"
#include "Framework/BaseTerm.hh"
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
template <typename CTERM, typename DTERM, int DIM>
class ArcJetLTEPhysicalModel : 
	public Framework::ConvectionDiffusionReactionPM
<CTERM, ArcJetTerm<DTERM>, ArcJetTerm<Framework::BaseTerm> > {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
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
    if (getNbEquations() == DIM+4) {
      CFLog(VERBOSE, "ArcJetLTEPhysicalModel::getTypeName() => ArcJetSALTE\n");
      return std::string("ArcJetSALTE" + Common::StringOps::to_str(DIM) + "D");
    } 
    cf_assert(getNbEquations() == DIM+3);
    CFLog(VERBOSE, "ArcJetLTEPhysicalModel::getTypeName() => ArcJetLTE\n");
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
  
private:
  
  /// number of turbulence equations
  CFuint m_nbTurbEqs;
  
}; // end of class ArcJetLTEPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetLTEPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_ArcJetLTEPhysicalModel_hh
