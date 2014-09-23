#ifndef COOLFluiD_Physics_StructMech_StructMech2DInertiaVarSet_hh
#define COOLFluiD_Physics_StructMech_StructMech2DInertiaVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/InertiaVarSet.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/NotImplementedException.hh"
#include "StructMechPhysicalModel.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the Inertia StructMech physical model 2D
   *
   * @author Thomas Wuilbaut
   */
class StructMech2DInertiaVarSet : public Framework::InertiaVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMech2D
   */
  StructMech2DInertiaVarSet(const std::string& name) :
    Framework::InertiaVarSet(name),
    _model(CFNULL)
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMech2DInertiaVarSet()
  {
  }

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    InertiaVarSet::configure(args);

    // add here configuration, specific of this class
  }

  /**
   * Setup the member data
   */
  virtual void setup()
  {
    CFAUTOTRACE;

    InertiaVarSet::setup();

    _model = Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<StructMechPhysicalModel>();
  }

  /**
   * Get the model
   */
   Common::SafePtr<StructMechPhysicalModel> getModel() const
   {
     return _model;
   }

private:

  /// physical model
  Common::SafePtr<StructMechPhysicalModel> _model;

}; // end of class StructMech2DInertiaVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech2DInertiaVarSet_hh
