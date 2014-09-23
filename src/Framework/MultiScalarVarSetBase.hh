// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MultiScalarVarSetBase_hh
#define COOLFluiD_Framework_MultiScalarVarSetBase_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/EquationSetData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BaseTerm;
    class State;
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a convective variable set with multi scalar equations
 * to append to an existing convective variable set
 *
 * @author Andrea Lani
 */
template <class BASEVS>
class MultiScalarVarSetBase : public BASEVS {
public: // classes
  typedef Framework::MultiScalarTerm<typename BASEVS::PTERM> PTERM;
  typedef typename BASEVS::EULERSET EULERSET;
  
  /**
   * Constructor
   */
  MultiScalarVarSetBase(Common::SafePtr<Framework::BaseTerm> term) :
    BASEVS(term),
    _mScalarModel(term.template d_castTo<PTERM>())
  {
  }
  
  /**
   * Default destructor
   */
  virtual ~MultiScalarVarSetBase()
  {
  }

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const = 0;

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
			     const RealVector& normal) 
  {
    BASEVS::splitJacobian(jacobPlus, jacobMin, eValues, normal); 
  }
  
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv,
					 RealVector& eValues,
					 const RealVector& normal) 
  {
    BASEVS::computeEigenValuesVectors(rightEv, leftEv, eValues, normal); 
  }
  
  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    BASEVS::setDimensionalValues(state, result);
  }
  
  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
				     RealVector& result)
  {
    BASEVS::setAdimensionalValues(state, result);
  }

  /// Set other dimensional values for useful physical quantities
  virtual void setDimensionalValuesPlusExtraValues
  (const State& state, RealVector& result, RealVector& extra)
  {
    BASEVS::setDimensionalValuesPlusExtraValues(state, result, extra);
  }
  
 /**
   * Get the model
   */
  Common::SafePtr<PTERM> getModel() const
  {
    cf_assert(_mScalarModel.isNotNull());
    return _mScalarModel;
  }

  /**
   * Get some data corresponding to the subset of equations related with
   * this variable set
   * @pre The most concrete ConvectiveVarSetBase will have to set these data
   */
  static std::vector<Framework::EquationSetData>& getEqSetData()
  {
    static std::vector<Framework::EquationSetData> eqSetData;
    return eqSetData;
  }
 
  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues) = 0;
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal) = 0;
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal) = 0;
  
protected:
  
  /**
   * Add an equation data with the specified number of equations
   */
  void addEqSetData(CFuint nbEqs)
  {
    using namespace std;
    
    vector<EquationSetData>& v = this->getEqSetData();
    const CFuint oldSize = v.size();
    CFLog(VERBOSE, "MultiScalarVarSetBase::addEqSetData() => \n");
    CFLog(VERBOSE, "MultiScalarVarSetBase::addEqSetData() => oldSize =" << oldSize << "\n");
    
    CFuint countEqs = EULERSET::getEqSetData()[0].size();
    const CFuint startSize = (countEqs > 0) ? oldSize + 1 : oldSize; 
    for (CFuint i = 0; i < oldSize; ++i) {
      countEqs += v[i].getEqSetVarIDs().size();
    }
    
    CFLog(VERBOSE, "MultiScalarVarSetBase::addEqSetData() => countEqs =" << countEqs << "\n");
    
    // set the equation set data for each of the equation subsets
    v.push_back(EquationSetData());
    v[oldSize].setup(startSize, countEqs, nbEqs);
  }
  
  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals) = 0;
  
  /// Computes the physical convective flux
  virtual void computeStateFlux(const RealVector& pdata) = 0;
  
  /**
   * Tells if the current equation sub set is considered
   */
  bool isCurrEqSS(CFuint iEqSS, CFuint size,
		  const std::vector<Framework::EquationSetData>& eqData)
  {
    for (CFuint i = 0; i < size; ++i) {
      if (iEqSS == eqData[i].getEqSetID()) return true;
    }
    return false;
  }
  
private:

  /// acquaintance of the model
  Common::SafePtr<PTERM> _mScalarModel;

}; // end of class MultiScalarVarSetBase

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MultiScalarVarSetBase_hh
