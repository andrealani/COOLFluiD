#ifndef COOLFluiD_Physics_StructMech_StructMech2DDisp_hh
#define COOLFluiD_Physics_StructMech_StructMech2DDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VectorialFunction.hh"
#include "StructMechPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a StructMech physical model 2D for conservative
 * variables
 *
 * @author Thomas Wuilbaut
 *
 */
class StructMech2DDisp : public Framework::ConvectiveVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMech2D
   */
  StructMech2DDisp(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~StructMech2DDisp();

  /**
   * Set up
   */
  void setup()
  {
    CFAUTOTRACE;
    //    Framework::ConvectiveVarSet::setup();

    _model = Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<StructMechPhysicalModel>();

  }

  /**
   * Get the physicalModel
   */
  Common::SafePtr<StructMechPhysicalModel> getModel() const
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
 
  /// Set the State correspoding to the PhysicalData
  virtual void computePhysicalData (const Framework::State& state, RealVector& pdata)
  {
    throw Common::NotImplementedException (FromHere(),"StructMech2DDisp::computePhysicalData()");
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
      (FromHere(), "StructMech2DDisp::computeFlux()");
  }

  /**
   * Compute the physical convective flux
   */
  void computeStateFlux(const RealVector& pdata)
  {
    throw Common::NotImplementedException
      (FromHere(), "StructMech2DDisp::computeStateFlux()");
  }

private:

  /// Corresponding physical model
  Common::SafePtr<StructMechPhysicalModel> _model;

}; // end of class StructMech2DDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech2DDisp_hh
