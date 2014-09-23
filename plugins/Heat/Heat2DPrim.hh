#ifndef COOLFluiD_Physics_Heat_Heat2DPrim_hh
#define COOLFluiD_Physics_Heat_Heat2DPrim_hh

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
 * This class represents a Heat physical model 2D for primitive
 * variables given in the following order { T }.
 *
 * @author Tiago Quintino
 *
 */
class Heat2DPrim : public Framework::ConvectiveVarSet {
public: // classes

  /**
   * Constructor
   * @see Heat2D
   */
  Heat2DPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Heat2DPrim();

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
  void computePhysicalData (const Framework::State& state, RealVector& pdata)
  {
   /// @todo broken after release 2009.3
   throw Common::NotImplementedException(FromHere(), "Heat2DPrim::computePhysicalData()");
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
  }

protected:

 /**
   * Compute the convective flux
   */
 void computeFlux(const RealVector& vars,
		   const RealVector& normals)
  {
    throw Common::NotImplementedException
      (FromHere(), "Heat2DPrim::computeFlux()");
  }

 /**
   * Compute the physical convective flux
   */
  void computeFlux(const Framework::State& vars)
  {
    throw Common::NotImplementedException
      (FromHere(), "Heat2DPrim::computeFlux()");
  }
  
  virtual void computeStateFlux(const RealVector& pdata)
  {
   throw Common::NotImplementedException
     (FromHere(), "Heat2DPrim::computeStateFlux()");
  }


private:

  /// Corresponding physical model
  Common::SafePtr<HeatPhysicalModel> _model;

}; // end of class Heat2DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat2DPrim_hh
