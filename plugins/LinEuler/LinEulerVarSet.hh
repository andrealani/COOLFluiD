/* WHAT THIS IS GOOD FOR?????????*/


#ifndef COOLFluiD_Physics_LinEuler_LinEulerVarSet_hh
#define COOLFluiD_Physics_LinEuler_LinEulerVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "LinEulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for Linearized Euler physical model.
 * @author Lilla Edit Koloszar
 * @author Nadege Villedieu
 * @author Tomas Kopacek
 * @author Matteo Parsani (modified for release 2009.3)
 */
class LinEulerVarSet : public Framework::ConvectiveVarSet {
public: // classes

  typedef LinEulerTerm PTERM;

  /**
   * Constructor
   * @see LinEulerPhysicalModel
   */
  LinEulerVarSet(Common::SafePtr<Framework::BaseTerm> term) :
    Framework::ConvectiveVarSet(term),
    _model(term.d_castTo<LinEulerTerm>()),
    _localMeanFlow(CFNULL)
  {
  }

  /**
   * Default destructor
   */
  virtual ~LinEulerVarSet()
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
    throw Common::NotImplementedException (FromHere(),"LinEulerVarSet::computeJacobians()");
  }

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
                          const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"LinEulerVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"LinEulerVarSet::computeEigenValuesVectors()");
  }

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
      (FromHere(), "LinEulerVarSet::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "LinEulerVarSet::setAdimensionalValues() not implemented");
  }

  /**
   * Get the model
   */
  Common::SafePtr<LinEulerTerm> getModel() const
  {
    return _model;
  }

  /**
   * Set the mean flow variables
   */
  // virtual void setExtraPhysicalVars(RealVector* extraVars)
  // {
  // cf_assert(extraVars->size() == getExtraPhysicalVarsSize());
  // _localMeanFlow = extraVars;
  // }

protected:

  /**
   * Set the vector of the eigenValues
   */
  virtual void computeEigenValues(const RealVector& pdata,
               const RealVector& normal, RealVector& result) = 0;


public:

  /// acquaintance of the model
  Common::SafePtr<LinEulerTerm> _model;

  /// pointer to a RealVector containing the local mean flow values
  RealVector* _localMeanFlow;

}; // end of class LinEulerVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_LinEulerVarSet_hh
