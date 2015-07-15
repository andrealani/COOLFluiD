#ifndef COOLFluiD_Physics_NavierStokes_EulerArcJetPhysicalModel_hh
#define COOLFluiD_Physics_NavierStokes_EulerArcJetPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/EulerTerm.hh"
#include "ArcJet/ArcJetTerm.hh"
#include "ArcJet/ArcJetInductionTerm.hh"
#include "Framework/ConvectionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an EulerArcJetPhysicalModel.
 *
 * @author Andrea Lani
 *
 */
template <int DIM>
class EulerArcJetPhysicalModel : public Framework::ConvectionReactionPM
<ArcJetInductionTerm<NavierStokes::EulerTerm>, 
 ArcJetTerm<Framework::BaseTerm> > {

public:
  
  /**
    * Constructor without arguments
    */
  EulerArcJetPhysicalModel(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~EulerArcJetPhysicalModel();
  
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
    return "Null";
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

}; // end of class EulerArcJetPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "EulerArcJetPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_EulerArcJetPhysicalModel_hh
