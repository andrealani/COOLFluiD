#ifndef COOLFluiD_Physics_Heat_Heat2DSourceVarSet_hh
#define COOLFluiD_Physics_Heat_Heat2DSourceVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SourceVarSet.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/NotImplementedException.hh"
#include "Heat/HeatPhysicalModel.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the Source Variable Set for Heat2D
   *
   * @author Thomas Wuilbaut
   */

class Heat2DSourceVarSet : public Framework::SourceVarSet {
public: // classes

  /**
   * Constructor
   * @see Heat2D
   */

  Heat2DSourceVarSet(const std::string& name) :
    Framework::SourceVarSet(name),
    _model(CFNULL)
  {
  }

  /**
   * Default destructor
   */
  virtual ~Heat2DSourceVarSet()
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
   * Setup the member data
   */
  virtual void setup()
  {
    Framework::SourceVarSet::setup();

    _model = Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<HeatPhysicalModel>();

  }

  /**
   * Gets the Linear Source Coeficients
   */
  virtual void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),"Heat2DSourceVarSet::getIndepSourceCoefs");
  }

  /**
   * Gets the Indep Source Coeficients
   */
  virtual void getIndepSourceCoefs(const Framework::State& state, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),"Heat2DSourceVarSet::getIndepSourceCoefs");
  }

  /**
   * Get the model
   */
   Common::SafePtr<HeatPhysicalModel> getModel() const
   {
     return _model;
   }


  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "Heat2DSourceVarSet";
  }


private:

  /// physical model
  Common::SafePtr<HeatPhysicalModel> _model;

}; // end of class Heat2DSourceVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat2DSourceVarSet_hh
