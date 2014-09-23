#ifndef COOLFluiD_Physics_ArcJet_ArcJetInductionDiffVarSet_hh
#define COOLFluiD_Physics_ArcJet_ArcJetInductionDiffVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "ArcJet/ArcJetInductionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework  {
   class PhysicalChemicalLibrary;
  }

  namespace Physics {
    
    namespace ArcJet {
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a ArcJet physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
template <typename BASEVS, typename ST>
class ArcJetInductionDiffVarSet : public BASEVS {
public: // classes

	typedef ArcJetInductionTerm<NavierStokes::EulerTerm> PTERM;

  /**
   * Constructor
   */
  ArcJetInductionDiffVarSet(const std::string& name,
                         Common::SafePtr<Framework::PhysicalModelImpl> model);


  /**
   * Default destructor
   */
  virtual ~ArcJetInductionDiffVarSet();
  
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
    throw Common::NotImplementedException (FromHere(),"ArcJetInductionDiffVarSet::getFlux()");
  }

  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);

private:
  
  /// reaction term
  Common::SafePtr<ST> m_arcJetReactionTerm;
  
  /// pointer to the physical chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

}; // end of class ArcJetInductionDiffVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetInductionDiffVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetInductionDiffVarSet_hh
