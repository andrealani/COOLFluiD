#ifndef COOLFluiD_Physics_Chemistry_CH4_ChemCH43DSourceVarSet_hh
#define COOLFluiD_Physics_Chemistry_CH4_ChemCH43DSourceVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SourceVarSet.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/NotImplementedException.hh"
#include "ChemCH4PhysicalModel.hh"
#include "Common/SafePtr.hh"

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
class ChemCH43DSourceVarSet : public Framework::SourceVarSet {
public:

  /**
   * Constructor
   */
  ChemCH43DSourceVarSet(const std::string& name) :
    Framework::SourceVarSet(name),
    _model(
      Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<ChemCH4PhysicalModel>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~ChemCH43DSourceVarSet()
  {
  }

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    Framework::SourceVarSet::configure(args);

    // add here configuration, specific of this class
  }

  /**
   * Gets the Linear Source Coeficients
   */
  virtual void getLinearSourceCoefs(const Framework::State& state, const RealVector& normals, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),"ChemCH43DSourceVarSet::getIndepSourceCoefs");
  }

  /**
   * Gets the Indep Source Coeficients
   */
  virtual void getIndepSourceCoefs(const Framework::State& state, const RealVector& normals, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),"ChemCH43DSourceVarSet::getIndepSourceCoefs");
  }

  /**
   * Get the model
   */
   Common::SafePtr<ChemCH4PhysicalModel> getModel() const
   {
     return _model;
   }


  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ChemCH43DSourceVarSet";
  }


private:

  /// physical model
  Common::SafePtr<ChemCH4PhysicalModel> _model;

}; // end of class ChemCH43DSourceVarSet

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_CH4_ChemCH43DSourceVarSet_hh
