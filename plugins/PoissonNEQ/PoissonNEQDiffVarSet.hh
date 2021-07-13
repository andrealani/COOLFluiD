#ifndef COOLFluiD_Physics_PoissonNEQ_PoissonNEQDiffVarSet_hh
#define COOLFluiD_Physics_PoissonNEQ_PoissonNEQDiffVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace PoissonNEQ {
      template <typename BASE> class PoissonNEQTerm;
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a PoissonNEQ physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
template <typename BASEVS>
class PoissonNEQDiffVarSet : public BASEVS {
public: // classes

  /**
   * Constructor
   */
  PoissonNEQDiffVarSet(const std::string& name,
                         Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~PoissonNEQDiffVarSet();

  /**
   * Set up private data
   */
  virtual void setup();
  
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
    throw Common::NotImplementedException (FromHere(),"PoissonNEQDiffVarSet::getFlux()");
  }

  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);

private:
  
  /// reaction term
  Common::SafePtr<PoissonNEQTerm<typename BASEVS::DTERM> > m_diffModel;

  /// electrical conductivity
  CFreal m_sigmaCoeff;

  /// array of tecmperatures
  RealVector m_tVec;
  
}; // end of class PoissonNEQDiffVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "PoissonNEQ/PoissonNEQDiffVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_PoissonNEQ_PoissonNEQDiffVarSet_hh
