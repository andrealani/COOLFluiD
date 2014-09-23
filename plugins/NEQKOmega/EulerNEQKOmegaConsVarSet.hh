#ifndef COOLFluiD_Physics_NEQKOmega_EulerNEQKOmegaConsVarSet_hh
#define COOLFluiD_Physics_NEQKOmega_EulerNEQKOmegaConsVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/EulerKOmegaConsVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a convective var set for k-Omega turbulence model
   *
   * @author Andrea Lani
   */
template <typename BASE, CFuint SGROUP>
class EulerNEQKOmegaConsVarSet : public KOmega::EulerKOmegaConsVarSet<BASE, SGROUP> {
public: // classes

  /**
   * Constructor
   */
  EulerNEQKOmegaConsVarSet(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~EulerNEQKOmegaConsVarSet();

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
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal);
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv,
					 RealVector& eValues,
					 const RealVector& normal);
  
  /**
   * Set other adimensional values for useful physical quantities
   */
  virtual void setDimensionalValuesPlusExtraValues
  (const Framework::State& state, RealVector& result, RealVector& extra);
  
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
  
  /// start ID for the turbulent equations
  CFuint m_startK;
  
  /// temporary array
  RealVector m_tmpResult;
  
  /// physical data
  RealVector m_pdatak;
  
}; // end of class EulerNEQKOmegaConsVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "EulerNEQKOmegaConsVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQKOmega_EulerNEQKOmegaConsVarSet_hh
