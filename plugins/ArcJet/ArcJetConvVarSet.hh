#ifndef COOLFluiD_Physics_ArcJet_ArcJetConvVarSet_hh
#define COOLFluiD_Physics_ArcJet_ArcJetConvVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    //   class BaseTerm;
    class State;
  }

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for LTE model and the electric potential
 *
 * @author Andrea Lani
 */
template <class BASEVS>
class ArcJetConvVarSet : public BASEVS {
public: // classes
  
  // typedef Framework::MultiScalarTerm<typename BASEVS::PTERM> PTERM;
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> PTERM;
  
  /**
   * Constructor
   */
  ArcJetConvVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~ArcJetConvVarSet();
 
  /// Set up the private data and give the maximum size of states physical
  /// data to store
  virtual void setup();
  
  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;

  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result);

  /**
   * Set other adimensional values for useful physical quantities
   */
  virtual void setDimensionalValuesPlusExtraValues(const Framework::State& state, 
						   RealVector& result,
						   RealVector& extra);
  
  /**
   * Get the model
   */
  Common::SafePtr<PTERM> getModel() const
  {
    cf_assert(_arcJetModel.isNotNull());
    return _arcJetModel;
  }
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  virtual void computePhysicalData(const Framework::State& state,
				   RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data,
					    Framework::State& state);
  
protected:

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const Framework::State& vars,
			   const RealVector& normals);


  /**
   * Compute the physical convective flux
   */
  virtual void computeFlux(const Framework::State& vars);
  
  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);
  
  virtual RealVector getFlux_array() const 
  {
    throw Common::NotImplementedException(FromHere(), "ArcJetConvVarSet::getFlux_array()");
  }
  
private:

  /// acquaintance of the model
  Common::SafePtr<PTERM> _arcJetModel;

}; // end of class ArcJetConvVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetConvVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetConvVarSet_hh
