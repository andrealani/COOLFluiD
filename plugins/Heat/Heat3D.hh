#ifndef COOLFluiD_Physics_Heat_Heat3D_hh
#define COOLFluiD_Physics_Heat_Heat3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "HeatPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a Heat3D.
 *
 * @author Tiago Quintino
 *
 */
class Heat3D : public HeatPhysicalModel {
public:

  /**
   * Constructor without arguments
   */
  Heat3D(const std::string& name);

  /**
   * Default destructor
   */
  ~Heat3D();
  
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

}; // end of class Heat3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat3D_hh
