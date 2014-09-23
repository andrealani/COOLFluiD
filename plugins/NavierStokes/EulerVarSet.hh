#ifndef COOLFluiD_Physics_NavierStokes_EulerVarSet_hh
#define COOLFluiD_Physics_NavierStokes_EulerVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/EquationSetData.hh"
#include "EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for an Euler physical model.
 *
 * @author Andrea Lani
 */
class EulerVarSet : public Framework::ConvectiveVarSet {
public: // classes
  typedef EulerVarSet EULERSET;
  typedef EulerTerm PTERM;
  
  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  EulerVarSet(Common::SafePtr<Framework::BaseTerm> term) :
    Framework::ConvectiveVarSet(term),
    _model(term.d_castTo<EulerTerm>()),
    _skipEnergyData(false),
    _physDataNeedCoordinates(false)
  {
  }
  
  /**
   * Default destructor
   */
  virtual ~EulerVarSet()
  {
  }

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();
  }

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const = 0;
   
  /**
   * Set the flag telling not to compute energy, enthalpy, sound speed in physical data 
   */
  void skipEnergyData(bool flag)
  {
    _skipEnergyData = flag;
  }
  
  /**
   * Set the jacobians
   */
  virtual void computeJacobians()
  {
    throw Common::NotImplementedException (FromHere(),"EulerVarSet::computeJacobians()");
  }
  
  /**
   * Compute the pressure derivative
   */
  virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp)
  {
    throw Common::NotImplementedException (FromHere(),"EulerVarSet::computePressureDerivatives()");
  }
  
  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
                          const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"EulerVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"EulerVarSet::computeEigenValuesVectors()");
  }
  
  /// Set the PhysicalData corresponding to the given State
  virtual void computePhysicalData (const Framework::State& state, RealVector& pdata) = 0;
  
  /**
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const = 0;

  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "EulerVarSet::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "EulerVarSet::setAdimensionalValues() not implemented");
  }

  /**
   * Get the model
   */
  Common::SafePtr<EulerTerm> getModel() const
  {
    return _model;
  }
  
  /**
   * Get some data corresponding to the subset of equations related with
   * this variable set
   * @pre The most concrete ConvectiveVarSet will have to set these data
   */
  static std::vector<Framework::EquationSetData>& getEqSetData()
  {
    static std::vector<Framework::EquationSetData> eqSetData;
    return eqSetData;
  }
  
private:
  
  /// acquaintance of the model
  Common::SafePtr<EulerTerm> _model;
  
protected:  
  
  /// flag telling not to compute energy, enthalpy, sound speed in physical data
  bool _skipEnergyData;
  
  /// flag telling if coordinates have to be stored in the physical data
  bool _physDataNeedCoordinates;
  
}; // end of class EulerVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_EulerVarSet_hh
