#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInletMHD2DProjection_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInletMHD2DProjection_hh

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
class SuperInletMHD2DProjection : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperInletMHD2DProjection(const std::string& name);

  /**
   * Default destructor
   */
  ~SuperInletMHD2DProjection();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private: //data

  /// phi value that is to be fixed
  CFreal _refPhi;

}; // end of class SuperInletMHD2DProjection

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInletMHD2DProjection_hh
