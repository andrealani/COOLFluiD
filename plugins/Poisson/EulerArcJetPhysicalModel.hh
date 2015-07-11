#ifndef COOLFluiD_Physics_Poisson_PoissonPhysicalModel_hh
#define COOLFluiD_Physics_Poisson_PoissonPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "ArcJet/PoissonTerm.hh"
#include "ArcJet/PoissonConvTerm.hh"
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
<PoissonConvTerm, PoissonTerm<Framework::BaseTerm> > {

public:
  
  /**
    * Constructor without arguments
    */
  PoissonPhysicalModel(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~PoissonPhysicalModel();
  
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

}; // end of class PoissonPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "PoissonPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonPhysicalModel_hh
