#ifndef COOLFluiD_Physics_Chemistry_CH4_ChemCH43DDiffusiveVarSet_hh
#define COOLFluiD_Physics_Chemistry_CH4_ChemCH43DDiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "Common/SafePtr.hh"
#include "ChemCH4PhysicalModel.hh"
#include "ChemCH43DPrim.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the diffusive Chemistry physical model 3D
 *
 * @author Tiago Quintino
 */
class ChemCH43DDiffusiveVarSet : public Framework::DiffusiveVarSet {
public:

  /**
   * Constructor
   */
  ChemCH43DDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    Framework::DiffusiveVarSet(name, model),
    _model(model.d_castTo<ChemCH4PhysicalModel>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~ChemCH43DDiffusiveVarSet()
  {
  }

  /**
   * Get the model
   */
   Common::SafePtr<ChemCH4PhysicalModel> getModel() const
   {
     return _model;
   }

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    Framework::DiffusiveVarSet::configure(args);

    // add here configuration, specific of this class
  }

private:

  /// physical model
  Common::SafePtr<ChemCH4PhysicalModel> _model;

}; // end of class ChemCH43DDiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_CH4_ChemCH43DDiffusiveVarSet_hh
