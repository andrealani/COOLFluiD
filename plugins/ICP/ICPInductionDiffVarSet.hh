#ifndef COOLFluiD_Physics_ICP_ICPInductionDiffVarSet_hh
#define COOLFluiD_Physics_ICP_ICPInductionDiffVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace ICP {
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a ICP physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
template <typename BASEVS, typename ST>
class ICPInductionDiffVarSet : public BASEVS {
public: // classes

  /**
   * Constructor
   */
  ICPInductionDiffVarSet(const std::string& name,
                         Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~ICPInductionDiffVarSet();
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius) 
  {
    throw Common::NotImplementedException (FromHere(),"ICPInductionDiffVarSet::getFlux()");
  }

  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);

private:
  
  /// reaction term
  Common::SafePtr<ST> m_icpReactionTerm;
  
}; // end of class ICPInductionDiffVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICPInductionDiffVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPInductionDiffVarSet_hh
