#ifndef COOLFluiD_Physics_StructMech_StructMech2DDiffusiveVarSet_hh
#define COOLFluiD_Physics_StructMech_StructMech2DDiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "StructMech2D.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the diffusive StructMech physical model 2D
   *
   * @author Thomas Wuilbaut
   */
class StructMech2DDiffusiveVarSet : public Framework::DiffusiveVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMech2D
   */

  StructMech2DDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model.d_castTo<StructMech2D>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMech2DDiffusiveVarSet()
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
   Common::SafePtr<StructMech2D> getModel() const
   {
     return _model;
   }


private:

  /// physical model
  Common::SafePtr<StructMech2D> _model;

}; // end of class StructMech2DDiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech2DDiffusiveVarSet_hh
