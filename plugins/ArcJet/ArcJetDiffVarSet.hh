#ifndef COOLFluiD_Physics_ArcJet_ArcJetDiffVarSet_hh
#define COOLFluiD_Physics_ArcJet_ArcJetDiffVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework  {
   class PhysicalChemicalLibrary;
  }

  namespace Physics {
    namespace NavierStokes {
      class EulerTerm;
    }
    
    namespace ArcJet {
      template <typename BASE> class ArcJetTerm;
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a ArcJet physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
template <typename BASEVS>
class ArcJetDiffVarSet : public BASEVS {
public: // classes

  /**
   * Constructor
   */
  ArcJetDiffVarSet(const std::string& name,
		   Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */

  virtual ~ArcJetDiffVarSet();
  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
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
    throw Common::NotImplementedException (FromHere(),"ArcJetDiffVarSet::getFlux()");
  }

  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);
  
  /// Compute the transport properties
  virtual void computeTransportProperties(const RealVector& state,
					  const std::vector<RealVector*>& gradients,
					  const RealVector& normal);
  
private:
  
  /// pointer to the physical chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// convective model
  Common::SafePtr<NavierStokes::EulerTerm> m_eulerModelLTE;
  
  /// diffusive model
  Common::SafePtr<ArcJetTerm<typename BASEVS::DTERM> > m_diffModel;
  
  /// frozen value of sigma
  CFreal m_sigmaCoeff;
  
}; // end of class ArcJetDiffVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetDiffVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetDiffVarSet_hh
