#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInletProjection_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInletProjection_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic inlet command in 2D
   * for projection scheme
   *
   * @author Radka Keslerova
   *
   */
class SuperInletProjection : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperInletProjection(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SuperInletProjection();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
private: //data
  
  /// phi value that is to be fixed
  CFreal _refPhi;
  CFint _inletCoronalBC;
  CFint _Phi_divB_zero;
  CFint _Phi_divB_extrapolated;
  
  /// array specifying the IDs for which a special treatment has to be applied
  std::vector<CFuint> _projectionIDs;
  
}; // end of class SuperInletProjection

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInletProjection_hh
