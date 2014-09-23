#ifndef COOLFluiD_Physics_ArcJet_ArcJetInductionConvVarSet_hh
#define COOLFluiD_Physics_ArcJet_ArcJetInductionConvVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Common/NotImplementedException.hh"
#include "NavierStokes/EulerTerm.hh"
#include "ArcJet/ArcJetInductionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Framework {
class State;
}

namespace Physics {

namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for chemical non-equilibrium models
 *
 * @author Andrea Lani
 */
template <class BASEVS>
class ArcJetInductionConvVarSet : public BASEVS {
public: // classes

	typedef ArcJetInductionTerm<NavierStokes::EulerTerm> PTERM;

	/**
	 * Constructor
	 */
	ArcJetInductionConvVarSet(Common::SafePtr<Framework::BaseTerm> term);

	/**
	 * Default destructor
	 */
	virtual ~ArcJetInductionConvVarSet();

	/**
	 * Set up the private data and give the maximum size of states physical
	 * data to store
	 */
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

	/**
	 * Split the jacobian
	 */
	virtual void splitJacobian(RealMatrix& jacobPlus,
			RealMatrix& jacobMin,
			RealVector& eValues,
			const RealVector& normal);

protected:

	/// Computes the convective flux projected on a normal
	virtual void computeFlux(const RealVector& pdata, const RealVector& normals);

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
		throw Common::NotImplementedException(FromHere(), "ArcJetInductionConvVarSet::getFlux_array()");
	}

private:

	/// acquaintance of the model
	Common::SafePtr<PTERM> _arcJetModel;

	/// right eigenvectors matrix
	RealMatrix _rmatEv;

	/// left eigenvectors matrix
	RealMatrix _lmatEv;

	/// array of positive eigenvalues
	RealVector _eValP;

	/// array of negative eigenvalues
	RealVector _eValM;

}; // end of class ArcJetInductionConvVarSet

//////////////////////////////////////////////////////////////////////////////

} // namespace ArcJet

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetInductionConvVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetInductionConvVarSet_hh
