#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionEpsDiffPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionEpsDiffPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD3DProjectionDiffPrim.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a MHD physical model 3D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class MHD3DProjectionEpsDiffPrim : public MHD3DProjectionDiffPrim {
public: // classes

  /**
   * Constructor
   * @see MHD3DProjectionDiffVarSet
   */
  MHD3DProjectionEpsDiffPrim(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~MHD3DProjectionEpsDiffPrim();
    
}; // end of class MHD3DProjectionEpsDiffPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionEpsDiffPrim_hh
