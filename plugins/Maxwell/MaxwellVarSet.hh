#ifndef COOLFluiD_Physics_Maxwell_MaxwellVarSet_hh
#define COOLFluiD_Physics_Maxwell_MaxwellVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Maxwell/ConvMaxwellTerm.hh"
#include "Framework/EquationSetData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for an Maxwell physical model.
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 */
class MaxwellVarSet : public Framework::ConvectiveVarSet {
public: // classes
  typedef MaxwellVarSet MAXWELLSET;
  typedef ConvMaxwellTerm PTERM;
  
  /**
   * Constructor
   * @see MaxwellPhysicalModel
   */
  MaxwellVarSet(Common::SafePtr<Framework::BaseTerm> term) :
    Framework::ConvectiveVarSet(term),
    _model(term.d_castTo<ConvMaxwellTerm>())
  {
  }
  
  /**
   * Default destructor
   */
  virtual ~MaxwellVarSet()
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
   * Set the jacobians
   */
  virtual void computeJacobians()
  {
    throw Common::NotImplementedException (FromHere(),"MaxwellVarSet::computeJacobians()");
  }
  
  /**
   * Compute the pressure derivative
   */
  virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp)
  {
    throw Common::NotImplementedException (FromHere(),"MaxwellVarSet::computePressureDerivatives()");
  }
  
  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
                          const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"MaxwellVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"MaxwellVarSet::computeEigenValuesVectors()");
  }
  
  /// Set the PhysicalData corresponding to the given State
  virtual void computePhysicalData (const Framework::State& state, RealVector& pdata) = 0;
  
  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "MaxwellVarSet::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "MaxwellVarSet::setAdimensionalValues() not implemented");
  }

  /**
   * Get the normal speed
   */
  CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const
  {
    return 1.;
  }

  /**
   * Get the model
   */
  Common::SafePtr<ConvMaxwellTerm> getModel() const
  {
    return _model;
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
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
  
protected:  
  
  /// @returns the number of equations of this VarSet
  CFuint getNbEqs() const { return 6; std::cout << "getNbEqs()\n"; }
  
private:
  
  /// acquaintance of the model
  Common::SafePtr<ConvMaxwellTerm> _model;
      
}; // end of class MaxwellVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_MaxwellVarSet_hh
