#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDisp_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VectorialFunction.hh"
#include "StructMechHeat2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a StructMechHeat physical model 2D for conservative
 * variables
 *
 * @author Thomas Wuilbaut
 *
 */
class StructMechHeat2DDisp : public Framework::ConvectiveVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMechHeat2D
   */
  StructMechHeat2DDisp(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~StructMechHeat2DDisp();

  /**
   * Set up
   */
  void setup()
  {
    CFAUTOTRACE;
    //    Framework::ConvectiveVarSet::setup();

    _model = Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<StructMechHeat2D>();

  }

  /**
   * Get the physicalModel
   */
  Common::SafePtr<StructMechHeat2D> getModel() const
  {
    return _model;
  }

  /**
   * Dimensional Values of state + extra values
   */
  void setDimensionalValuesPlusExtraValues(const Framework::State& state,
    RealVector& result,RealVector& extra);

  /**
   * Get Names of Extra Variables
   */
  std::vector<std::string> getExtraVarNames() const;
  
  /// Set the PhysicalData corresponding to the given State
  void computePhysicalData (const Framework::State& state, RealVector& pdata);
 
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
  }
  
protected:

  /**
   * Compute the convective flux
   */
  void computeFlux(const Framework::State& vars,
		   const RealVector& normals)
  {
    throw Common::NotImplementedException
      (FromHere(), "StructMechHeat2DDisp::computeFlux()");
  }

  /**
   * Compute the physical convective flux
   */
  void computeFlux(const Framework::State& vars)
  {
    throw Common::NotImplementedException
      (FromHere(), "StructMechHeat2DDisp::computeFlux()");
  }
  
  /// @todo broken after release 2009.3
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  /// @todo broken after release 2009.3
  virtual void computeStateFlux(const RealVector& pdata);

private:

  /// Corresponding physical model
  Common::SafePtr<StructMechHeat2D> _model;

}; // end of class StructMechHeat2DDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDisp_hh
