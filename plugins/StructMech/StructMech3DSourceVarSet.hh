#ifndef COOLFluiD_Physics_StructMech_StructMech3DSourceVarSet_hh
#define COOLFluiD_Physics_StructMech_StructMech3DSourceVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalModel.hh"
#include "Framework/SourceVarSet.hh"
#include "Common/NotImplementedException.hh"
#include "StructMechPhysicalModel.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Source StructMech physical model 3D
 *
 * @author Tiago Quintino
 */
class StructMech3DSourceVarSet : public Framework::SourceVarSet {
public: // classes

  /**
   * Constructor
   */
  StructMech3DSourceVarSet(const std::string& name) :
    Framework::SourceVarSet(name),
    _model(
      Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<StructMechPhysicalModel>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~StructMech3DSourceVarSet()
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
  virtual void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),"Heat3DSourceVarSet::getIndepSourceCoefs");
  }

  /**
   * Gets the Linear Source Coeficients
   */
  virtual void getIndepSourceCoefs(const Framework::State& state, const RealVector& normals, RealMatrix& coef)
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

}; // end of class StructMech3DSourceVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech3DSourceVarSet_hh
