#ifndef COOLFluiD_Physics_ICPNEQ_ICPNEQPhysicalModel_hh
#define COOLFluiD_Physics_ICPNEQ_ICPNEQPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarTerm.hh"
#include "NEQ/NavierStokesNEQ.hh"
#include "NEQ/NEQReactionTerm.hh"
#include "ICP/ICPReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an ICPNEQPhysicalModel.
 *
 * @author Andrea Lani
 */
template <int DIM>
class ICPNEQPhysicalModel : public NEQ::NavierStokesNEQ
<DIM, 
 Framework::MultiScalarTerm<NavierStokes::EulerTerm>, 
 NavierStokes::NSTerm, 
 ICPReactionTerm<NEQ::NEQReactionTerm> > {

public:
  
  /**
   * Constructor without arguments
   */
  ICPNEQPhysicalModel(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ICPNEQPhysicalModel();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  virtual std::string getTypeName() const
  {
    return std::string("ICPNEQ" + Common::StringOps::to_str(DIM) + "D");
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

}; // end of class ICPNEQPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICPNEQPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICPNEQ_ICPNEQPhysicalModel_hh
