#ifndef COOLFluiD_Physics_StructMech_StructMech2DInertiaDisp_hh
#define COOLFluiD_Physics_StructMech_StructMech2DInertiaDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMech2DInertiaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a StructMech physical model 2D for Dispitive
   * variables
   *
   * @author Thomas Wuilbaut
   */

class StructMech2DInertiaDisp : public StructMech2DInertiaVarSet {
public:

  /**
   * Constructor
   */
  StructMech2DInertiaDisp(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMech2DInertiaDisp();

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

}; // end of class StructMech2DInertiaDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech2DInertia_hh
