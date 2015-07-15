#ifndef COOLFluiD_Physics_NavierStokes_ArcJetPhysicalModel_hh
#define COOLFluiD_Physics_NavierStokes_ArcJetPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/NSTerm.hh"
#include "ArcJet/ArcJetTerm.hh"
#include "ArcJet/ArcJetInductionTerm.hh"
#include "Framework/ConvectionDiffusionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an ArcJetPhysicalModel.
 *
 * @author Andrea Lani
 *
 */
template <int DIM>
class ArcJetPhysicalModel : public Framework::ConvectionDiffusionReactionPM<
  ArcJetInductionTerm<NavierStokes::EulerTerm>, NavierStokes::NSTerm, 
  ArcJetTerm<Framework::BaseTerm> > {

public:
  
  /**
    * Constructor without arguments
    */
  ArcJetPhysicalModel(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ArcJetPhysicalModel();
  
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
    return std::string("ArcJet" + Common::StringOps::to_str(DIM) + "D");
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

}; // end of class ArcJetPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_ArcJetPhysicalModel_hh
