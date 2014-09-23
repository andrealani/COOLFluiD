#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat2DInertiaDisp_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat2DInertiaDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMechHeat2DInertiaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a StructMechHeat physical model 2D for Dispitive
   * variables
   *
   * @author Thomas Wuilbaut
   */

class StructMechHeat2DInertiaDisp : public StructMechHeat2DInertiaVarSet {
public:

  /**
   * Constructor
   */
  StructMechHeat2DInertiaDisp(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMechHeat2DInertiaDisp();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Setup the member data
   */
  virtual void setup();

private:

  /// Storage of the Density
  CFreal _rho;

}; // end of class StructMechHeat2DInertiaDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeat2DInertia_hh
