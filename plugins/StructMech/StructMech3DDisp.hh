#ifndef COOLFluiD_Physics_StructMech_StructMech3DDisp_hh
#define COOLFluiD_Physics_StructMech_StructMech3DDisp_hh

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
 * This class represents a StructMech physical model 3D for conservative
 * variables
 *
 * @author Tiago Quintino
 *
 */
class StructMech3DDisp : public Framework::ConvectiveVarSet {
public: // classes

  /**
   * Constructor
   * @see StructMech3D
   */
  StructMech3DDisp(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~StructMech3DDisp();

  /**
   * Set up
   */
  void setup();

  /**
   * Get the physicalModel
   */
  Common::SafePtr<StructMechPhysicalModel> getModel() const
  {
    return _model;
  }
  
  /// Set the State correspoding to the PhysicalData
  virtual void computePhysicalData (const Framework::State& state, RealVector& pdata)
  {
    throw Common::NotImplementedException (FromHere(),"StructMech3DDisp::computePhysicalData()");
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
      (FromHere(), "StructMech3DDisp::computeFlux()");
  }

  /**
   * Compute the physical convective flux
   */
  void computeStateFlux(const RealVector& pdata)
  {
    throw Common::NotImplementedException
      (FromHere(), "StructMech3DDisp::computeStateFlux()");
  }

private:

  /// Corresponding physical model
  Common::SafePtr<StructMechPhysicalModel> _model;

}; // end of class StructMech3DDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech3DDisp_hh
