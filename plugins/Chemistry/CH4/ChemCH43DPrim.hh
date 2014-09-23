#ifndef COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrim_hh
#define COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VectorialFunction.hh"
#include "ChemCH4PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

enum PrimIndex { XCH4, XO2, XCO2, XH2O };

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a ChemCH4 physical model 3D for primitive
 * variables given in the following order { XCH4, XO2, XCO2, XCO }.
 *
 * @author Tiago Quintino
 *
 */
class ChemCH43DPrim : public Framework::ConvectiveVarSet {
public: // classes

  /**
   * Constructor
   * @see ChemCH43D
   */
  ChemCH43DPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~ChemCH43DPrim();

  /**
   * Set up
   */
  void setup()
  {
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs) 
  {
  }
  
  /**
   * Get the physicalModel
   */
  Common::SafePtr<ChemCH4PhysicalModel> getModel() const
  {
    return _model;
  }
  
  
private: // methods

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& vars,
                          const RealVector& normals);

  /**
   * Compute the physical convective flux
   */
  virtual void computeFlux(const RealVector& vars);
  
  virtual void computeStateFlux(const RealVector& vars);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  
  void computePhysicalData(const Framework::State& state, RealVector& data);

 
  /**
   * Set the PhysicalData corresponding to the given State
   */
  void computeStateFromPhysicalData(const RealVector& data, Framework::State& state)
  {
  }

  /**
   * Set the vector of the eigenValues
   */
  void computeEigenValues_Impl(Framework::State& state,
			  const RealVector& normal, RealVector& eValues)
  {
  }

  /**
   * Get the maximum eigen value
   */
  CFreal getMaxEigenValueImpl(Framework::State& state, const RealVector& normal)
  {
    return 0.0;
  }

  /**
   * Get the maximum absolute eigenvalue
   */
  CFreal getMaxAbsEigenValueImpl(Framework::State& state, const RealVector& normal)
  {
    return 0.0;
  }

private: // data

  /// Corresponding physical model
  Common::SafePtr<ChemCH4PhysicalModel> _model;

}; // end of class ChemCH43DPrim

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_ChemCH43DPrim_hh
