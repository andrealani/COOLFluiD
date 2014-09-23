#ifndef COOLFluiD_Numerics_FiniteVolume_SteadySilverMullerCondition2D_hh
#define COOLFluiD_Numerics_FiniteVolume_SteadySilverMullerCondition2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic inlet command
   *
   * @author Thomas Wuilbaut
   *
   */
class SteadySilverMullerCondition2D : public SuperInlet {
public:

  /**
   * Constructor
   */
  SteadySilverMullerCondition2D(const std::string& name);

  /**
   * Default destructor
   */
  ~SteadySilverMullerCondition2D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

protected: // data

  /// storage for the temporary boundary point coordinates
  RealVector _variables;

}; // end of class SteadySilverMullerCondition2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SteadySilverMullerCondition2D_hh
