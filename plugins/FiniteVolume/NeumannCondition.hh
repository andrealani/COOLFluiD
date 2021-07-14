#ifndef COOLFluiD_Numerics_FiniteVolume_NeumannCondition_hh
#define COOLFluiD_Numerics_FiniteVolume_NeumannCondition_hh

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
   * @author Alejandro Alvarez
   *
   */
class NeumannCondition : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  // static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NeumannCondition(const std::string& name);

  /**
   * Default destructor
   */
  ~NeumannCondition();

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
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

protected: // data

  /// storage for the temporary boundary point coordinates
  RealVector _variables;
  
/// direction vector for RL
  RealVector _eRL;
  
  /// value at the boundary
  CFreal _value;

}; // end of class NeumannCondition

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NeumannCondition_hh
