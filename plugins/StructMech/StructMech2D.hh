#ifndef COOLFluiD_Physics_StructMech_StructMech2D_hh
#define COOLFluiD_Physics_StructMech_StructMech2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMechPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a StructMech2D.
 *
 * @author Thomas Wuilbaut
 *
 */
class StructMech2D : public StructMechPhysicalModel {
public:

  /**
   * Constructor without arguments
   */
  StructMech2D(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMech2D();

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * Check if this state is in a valid state
   */
  bool validate(const Framework::State& state) const;

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

private:

}; // end of class StructMech2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech2D_hh
