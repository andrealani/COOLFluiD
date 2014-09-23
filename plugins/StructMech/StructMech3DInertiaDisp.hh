#ifndef COOLFluiD_Physics_StructMech_StructMech3DInertiaDisp_hh
#define COOLFluiD_Physics_StructMech_StructMech3DInertiaDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMech3DInertiaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a StructMech physical model 3D for Dispitive
   * variables
   *
   * @author Thomas Wuilbaut
   */
class StructMech3DInertiaDisp : public StructMech3DInertiaVarSet {
public:

  /**
   * Constructor
   */
  StructMech3DInertiaDisp(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMech3DInertiaDisp();

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

}; // end of class StructMech3DInertiaDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech3DInertia_hh
