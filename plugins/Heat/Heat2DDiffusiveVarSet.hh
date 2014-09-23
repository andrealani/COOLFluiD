#ifndef COOLFluiD_Physics_Heat_Heat2DDiffusiveVarSet_hh
#define COOLFluiD_Physics_Heat_Heat2DDiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "Heat2D.hh"
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

class Heat2DDiffusiveVarSet : public Framework::DiffusiveVarSet {
public:

  /**
   * Constructor
   */
  Heat2DDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model.d_castTo<Heat2D>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~Heat2DDiffusiveVarSet()
  {
  }

  /**
   * Get the model
   */
   Common::SafePtr<Heat2D> getModel() const
   {
     return _model;
   }

private:

  /// physical model
  Common::SafePtr<Heat2D> _model;

}; // end of class Heat2DDiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat2DDiffusiveVarSet_hh
