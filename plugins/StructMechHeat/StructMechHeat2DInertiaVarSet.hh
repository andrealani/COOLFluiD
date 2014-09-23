#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat2DInertiaVarSet_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat2DInertiaVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/InertiaVarSet.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/NotImplementedException.hh"
#include "StructMechHeatPhysicalModel.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the Inertia StructMechHeat physical model 2D
   *
   * @author Thomas Wuilbaut
   */
class StructMechHeat2DInertiaVarSet : public Framework::InertiaVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMechHeat2D
   */
  StructMechHeat2DInertiaVarSet(const std::string& name) :
    Framework::InertiaVarSet(name),
    _model(CFNULL)
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMechHeat2DInertiaVarSet()
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
      getImplementor().d_castTo<StructMechHeatPhysicalModel>();
  }

  /**
   * Get the model
   */
   Common::SafePtr<StructMechHeatPhysicalModel> getModel() const
   {
     return _model;
   }

private:

  /// physical model
  Common::SafePtr<StructMechHeatPhysicalModel> _model;

}; // end of class StructMechHeat2DInertiaVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeat2DInertiaVarSet_hh
