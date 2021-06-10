#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionDiffVarSet_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionDiffVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a MHD physical model 3D for generic variables
   *
   * @author Andrea Lani
   */
class MHD3DProjectionDiffVarSet : public MHDProjectionDiffVarSet {
public: // classes

  /**
   * Constructor
   * @see MHD3D
   */
  MHD3DProjectionDiffVarSet(const std::string& name,
			    Common::SafePtr<Framework::PhysicalModelImpl> model) :
    MHDProjectionDiffVarSet(name, model)
  {
  }

  /**
   * Default destructor
   */
  virtual ~MHD3DProjectionDiffVarSet()
  {
  }

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
  virtual RealMatrix& getFlux(const RealVector& state,
			      const std::vector<RealVector*>& gradients,
			      const CFreal& radius);
  
}; // end of class MHD3DProjectionDiffVarSet
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionDiffVarSet_hh
