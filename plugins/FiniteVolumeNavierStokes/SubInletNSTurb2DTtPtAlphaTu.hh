#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb2DTtPtAlphaTu_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb2DTtPtAlphaTu_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "SubInletEuler2DTtPtAlpha.hh"
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with the initial conditions given for tTotal, pTotal and alpha
   *
   * @author Khalil Bensassi 
   *
   */
class SubInletNSTurb2DTtPtAlphaTu : public SubInletEuler2DTtPtAlpha  {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletNSTurb2DTtPtAlphaTu(const std::string& name);

  /**
   * Default destructor
   */
  virtual  ~SubInletNSTurb2DTtPtAlphaTu();

    /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

 
  protected: 
  CFreal     _Tu;

}; // end of class SubInletNSTurb2DTtPtAlphaTu

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb2DTtPtAlphaTu_hh
