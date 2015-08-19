#ifndef COOLFluiD_Numerics_FiniteVolume_DirichletCondition_hh
#define COOLFluiD_Numerics_FiniteVolume_DirichletCondition_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Dirichlet condition command
   *
   * @author Alejandro Alvarez
   *
   */
class DirichletCondition : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  DirichletCondition(const std::string& name);

  /**
   * Default destructor
   */
  ~DirichletCondition();

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

  /// a string holding the function definition
  std::vector<std::string> _function;

  /// storage for the temporary boundary point coordinates
  RealVector _variables;

  /// value at the boundary
  CFreal _value;

}; // end of class DirichletCondition

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DirichletCondition_hh
