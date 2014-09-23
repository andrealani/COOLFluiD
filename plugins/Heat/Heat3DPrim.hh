#ifndef COOLFluiD_Physics_Heat_Heat3DPrim_hh
#define COOLFluiD_Physics_Heat_Heat3DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VectorialFunction.hh"
#include "HeatPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Heat physical model 3D for primitive
 * variables given in the following order { T }.
 *
 * @author Tiago Quintino
 *
 */
class Heat3DPrim : public Framework::ConvectiveVarSet {
public: // classes

  /**
   * Constructor
   * @see Heat3D
   */
  Heat3DPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Heat3DPrim();

  /**
   * Set up
   */
  void setup()
  {
    _model = Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<HeatPhysicalModel>();
  }

  /**
   * Get the physicalModel
   */
  Common::SafePtr<HeatPhysicalModel> getModel() const
  {
    return _model;
  }

  /// Set the PhysicalData corresponding to the given State
  void computePhysicalData(const Framework::State & state, RealVector& pdata)
  {
   throw Common::NotImplementedException
     (FromHere(), "Heat3DPrim::computePhysicalData()");
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3;
  }
  
protected:

  /**
   * Compute the convective flux
   */
  void computeFlux(const RealVector& vars,
		   const RealVector& normals)
  {
    throw Common::NotImplementedException
      (FromHere(), "Heat3DPrim::computeFlux()");
  }
  
  void computeStateFlux(const RealVector& vars)
  {
   /// @todo broken after release 2009.3
   throw Common::NotImplementedException
     (FromHere(), "Heat3DPrim::computeStateFlux()");
  }

  /**
   * Compute the phys convective flux
   */
  void computeFlux(const Framework::State& vars)
  {
    throw Common::NotImplementedException
      (FromHere(), "Heat3DPrim::computeFlux()");
  }

private:

  /// Corresponding physical model
  Common::SafePtr<HeatPhysicalModel> _model;

}; // end of class Heat3DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat3DPrim_hh
