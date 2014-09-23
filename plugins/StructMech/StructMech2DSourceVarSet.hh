#ifndef COOLFluiD_Physics_StructMech_StructMech2DSourceVarSet_hh
#define COOLFluiD_Physics_StructMech_StructMech2DSourceVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SourceVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "StructMechPhysicalModel.hh"
#include "Common/SafePtr.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the Source StructMech physical model 2D
   *
   * @author Thomas Wuilbaut
   */
class StructMech2DSourceVarSet : public Framework::SourceVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMech2D
   */

  StructMech2DSourceVarSet(const std::string& name) :
    Framework::SourceVarSet(name),
    _model(CFNULL)
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMech2DSourceVarSet()
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
      getImplementor().d_castTo<StructMechPhysicalModel>();

  }

  /**
   * Gets the Linear Source Coeficients
   */
  virtual void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),"Heat3DSourceVarSet::getIndepSourceCoefs");
  }

  /**
   * Gets the Linear Source Coeficients
   */
  virtual void getIndepSourceCoefs(const Framework::State& state, const RealVector& normals, RealVector& coef)
  {
    throw Common::NotImplementedException (FromHere(),"Heat3DSourceVarSet::getIndepSourceCoefs");
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

}; // end of class StructMech2DSourceVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech2DSourceVarSet_hh
