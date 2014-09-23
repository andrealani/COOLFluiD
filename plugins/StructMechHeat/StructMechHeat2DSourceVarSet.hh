#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat2DSourceVarSet_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat2DSourceVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SourceVarSet.hh"
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
   * This class represents the Source StructMechHeat physical model 2D
   *
   * @author Thomas Wuilbaut
   */
class StructMechHeat2DSourceVarSet : public Framework::SourceVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMechHeat2D
   */

  StructMechHeat2DSourceVarSet(const std::string& name) :
    Framework::SourceVarSet(name),
    _model(CFNULL)
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMechHeat2DSourceVarSet()
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
    CFAUTOTRACE;

    Framework::SourceVarSet::setup();

    _model = Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<StructMechHeatPhysicalModel>();

  }

  /**
   * Gets the Linear Source Coeficients
   */
  virtual void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),"StructMechHeat2DSourceVarSet::getIndepSourceCoefs");
  }

  /**
   * Gets the Linear Source Coeficients
   */
  virtual void getIndepSourceCoefs(const Framework::State& state, const RealVector& normals, RealVector& coef)
  {
    throw Common::NotImplementedException (FromHere(),"StructMechHeat2DSourceVarSet::getIndepSourceCoefs");
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

}; // end of class StructMechHeat2DSourceVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeat2DSourceVarSet_hh
