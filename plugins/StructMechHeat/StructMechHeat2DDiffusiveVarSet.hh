#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDiffusiveVarSet_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "StructMechHeat2D.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the diffusive StructMechHeat physical model 2D
   *
   * @author Thomas Wuilbaut
   */
class StructMechHeat2DDiffusiveVarSet : public Framework::DiffusiveVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMechHeat2D
   */

  StructMechHeat2DDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model.d_castTo<StructMechHeat2D>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMechHeat2DDiffusiveVarSet()
  {
  }

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    CFAUTOTRACE;

    DiffusiveVarSet::configure(args);

    // add here configuration, specific of this class
  }

  /**
   * Get the model
   */
   Common::SafePtr<StructMechHeat2D> getModel() const
   {
     return _model;
   }


private:

  /// physical model
  Common::SafePtr<StructMechHeat2D> _model;

}; // end of class StructMechHeat2DDiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDiffusiveVarSet_hh
