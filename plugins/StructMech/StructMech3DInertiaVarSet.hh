#ifndef COOLFluiD_Physics_StructMech_StructMech3DInertiaVarSet_hh
#define COOLFluiD_Physics_StructMech_StructMech3DInertiaVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/InertiaVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/PhysicalModel.hh"
#include "StructMechPhysicalModel.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the Inertia StructMech physical model 3D
   *
   * @author Thomas Wuilbaut
   */
class StructMech3DInertiaVarSet : public Framework::InertiaVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMech3D
   */

  StructMech3DInertiaVarSet(const std::string& name) :
    Framework::InertiaVarSet(name),
    _model(CFNULL)
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMech3DInertiaVarSet()
  {
  }

  /**
   * Setup the member data
   */
  virtual void setup()
  {
    InertiaVarSet::setup();

    _model = Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<StructMechPhysicalModel>();
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
   * Gets the Mass Matrix
   */
  virtual void getWeakMassMat(const CFreal& W,
                              const CFreal& N,
                              const Framework::State& state,
                              RealMatrix& result)
  {
    throw Common::NotImplementedException (FromHere(),"StructMech3DInertiaVarSet::getWeakMassMat()");
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

}; // end of class StructMech3DInertiaVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech3DInertiaVarSet_hh
