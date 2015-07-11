#ifndef COOLFluiD_Physics_Poisson_PoissonDiffVarSet_hh
#define COOLFluiD_Physics_Poisson_PoissonDiffVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Poisson/PoissonDiffTerm.hh"
#include "Poisson/PoissonConvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework  {
   class PhysicalChemicalLibrary;
  }

  namespace Physics {
    
    namespace Poisson {
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Poisson physical model 2D for primitive
   * variables
   *
   * @author Alejandro Alvarez
   */

class PoissonDiffVarSet : public Framework::DiffusiveVarSet {
public: // classes

  /**
   * Constructor
   */
  PoissonDiffVarSet(const std::string& name,
		   Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */

  virtual ~PoissonDiffVarSet();
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
    throw Common::NotImplementedException (FromHere(),"PoissonDiffVarSet::getFlux()");
  }

  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);

private:

  /// pointer to the physical chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
   
}; // end of class PoissonDiffVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonDiffVarSet_hh
