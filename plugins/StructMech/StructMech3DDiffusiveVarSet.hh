#ifndef COOLFluiD_Physics_StructMech_StructMech3DDiffusiveVarSet_hh
#define COOLFluiD_Physics_StructMech_StructMech3DDiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "StructMechPhysicalModel.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the diffusive StructMech physical model 3D
 *
 * @author Tiago Quintino
 */
class StructMech3DDiffusiveVarSet : public Framework::DiffusiveVarSet {
public: // classes

  /**
   * Constructor
   */
  StructMech3DDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model.d_castTo<StructMechPhysicalModel>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMech3DDiffusiveVarSet()
  {
  }

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    DiffusiveVarSet::configure(args);

    // add here configuration, specific of this class
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

}; // end of class StructMech3DDiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech3DDiffusiveVarSet_hh
