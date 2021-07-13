#ifndef COOLFluiD_Physics_PoissonNEQ_PoissonNEQPhysicalModel_hh
#define COOLFluiD_Physics_PoissonNEQ_PoissonNEQPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarTerm.hh"
#include "NEQ/NavierStokesNEQ.hh"
#include "NEQ/NEQReactionTerm.hh"
#include "PoissonNEQ/PoissonNEQTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an PoissonNEQPhysicalModel.
 *
 * @author Andrea Lani
 */
template <int DIM>
class PoissonNEQPhysicalModel : public NEQ::NavierStokesNEQ
<DIM, 
 Framework::MultiScalarTerm<NavierStokes::EulerTerm>, 
 PoissonNEQTerm<NavierStokes::NSTerm>, 
 NEQ::NEQReactionTerm> {
  
public:
  
  /**
   * Constructor without arguments
   */
  PoissonNEQPhysicalModel(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~PoissonNEQPhysicalModel();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  virtual std::string getTypeName() const
  {
    return std::string("PoissonNEQ" + Common::StringOps::to_str(DIM) + "D");
  }

  /**
   * Get the convective name
   */
  virtual std::string getConvectiveName() const;

  /**
   * Get the diffusive name
   */
  virtual std::string getDiffusiveName() const
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

}; // end of class PoissonNEQPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "PoissonNEQPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_PoissonNEQ_PoissonNEQPhysicalModel_hh
