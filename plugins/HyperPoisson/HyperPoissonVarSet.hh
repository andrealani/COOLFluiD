#ifndef COOLFluiD_Physics_HyperPoisson_HyperPoissonVarSet_hh
#define COOLFluiD_Physics_HyperPoisson_HyperPoissonVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/EquationSetData.hh"
#include "HyperPTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for an Hyperbolized Poisson physical model.
 *
 * @author Rayan Dhib
 * @author Andrea Lani
 */
class HyperPoissonVarSet : public Framework::ConvectiveVarSet {
public: // classes
  typedef HyperPTerm PTERM;
  
  /**
   * Constructor
   * @see HyperPoissonPhysicalModel
   */
  HyperPoissonVarSet(Common::SafePtr<Framework::BaseTerm> term) :
    Framework::ConvectiveVarSet(term),
    _model(term.d_castTo<HyperPTerm>()),
    _physDataNeedCoordinates(false)
  {
  }
  
  /**
   * Default destructor
   */
  virtual ~HyperPoissonVarSet()
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
    throw Common::NotImplementedException (FromHere(),"HyperPoissonVarSet::computeJacobians()");
  }
  
  /**
   * Compute the pressure derivative
   */
  virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp)
  {
    throw Common::NotImplementedException (FromHere(),"HyperPoissonVarSet::computePressureDerivatives()");
  }
  
  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
                          const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"HyperPoissonVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"HyperPoissonVarSet::computeEigenValuesVectors()");
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
      (FromHere(), "HyperPoissonVarSet::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "HyperPoissonVarSet::setAdimensionalValues() not implemented");
  }

  /**
   * Get the model
   */
  Common::SafePtr<HyperPTerm> getModel() const
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
  Common::SafePtr<HyperPTerm> _model;
  
protected:  
  
  /// flag telling if coordinates have to be stored in the physical data
  bool _physDataNeedCoordinates;
  
}; // end of class HyperPoissonVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_HyperPoisson_HyperPoissonVarSet_hh
