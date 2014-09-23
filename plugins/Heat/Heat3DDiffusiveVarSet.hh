#ifndef COOLFluiD_Physics_Heat_Heat3DDiffusiveVarSet_hh
#define COOLFluiD_Physics_Heat_Heat3DDiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "HeatPhysicalModel.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the diffusive Heat physical model 3D
   *
   * @author Thomas Wuilbaut
   */
class Heat3DDiffusiveVarSet : public Framework::DiffusiveVarSet {
public:

  /**
   * Constructor
   */
  Heat3DDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model.d_castTo<HeatPhysicalModel>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~Heat3DDiffusiveVarSet()
  {
  }

  /**
   * Get the model
   */
   Common::SafePtr<HeatPhysicalModel> getModel() const
   {
     return _model;
   }

private:

  /// physical model
  Common::SafePtr<HeatPhysicalModel> _model;

}; // end of class Heat3DDiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat3DDiffusiveVarSet_hh
