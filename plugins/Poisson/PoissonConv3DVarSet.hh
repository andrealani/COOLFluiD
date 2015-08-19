#ifndef COOLFluiD_Physics_Poisson_PoissonConv3DVarSet_hh
#define COOLFluiD_Physics_Poisson_PoissonConv3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "PoissonConvVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 3D Poisson physical model.
 *
 * @author Andrea Lani
 */
class PoissonConv3DVarSet : public PoissonConvVarSet {
public: // classes
  
  typedef PoissonConv3DVarSet EULERSET;
  
  /**
   * Constructor
   * @see PoissonPhysicalModel
   */
  PoissonConv3DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~PoissonConv3DVarSet();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

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
			  const RealVector& normal) = 0;

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal) = 0;
  
  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);
  
protected:

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& pdata,
			   const RealVector& normals);

  /**
   * Compute the convective flux
   */
  virtual void computeStateFlux(const RealVector& vars);
  
}; // end of class PoissonConv3DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonConv3DVarSet_hh
