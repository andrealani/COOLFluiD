#ifndef COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmegaVarSet_hh
#define COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmegaVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model for chemical NEQ
   * and KOmega rubulence models
   *
   * @author Andrea Lani
   */
template <typename BASE>
class NavierStokesNEQKOmegaVarSet : public BASE {
public: // classes
  
  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesNEQKOmegaVarSet(const std::string& name,
			Common::SafePtr<Framework::PhysicalModelImpl> model);


  /**
   * Default destructor
   */
  
  virtual ~NavierStokesNEQKOmegaVarSet();

  /**
   * Set up private data
   */
  virtual void setup();
  
  /// Compute the transport properties
  virtual void computeTransportProperties(const RealVector& state,
					  const std::vector<RealVector*>& gradients,
					  const RealVector& normal);
  
}; // end of class NavierStokesNEQKOmegaVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQKOmegaVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmegaVarSet_hh
