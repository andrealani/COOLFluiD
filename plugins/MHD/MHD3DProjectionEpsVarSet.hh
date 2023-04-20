#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionEpsVarSet_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionEpsVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD3DProjectionVarSet.hh"
#include "MHD/MHDProjectionEpsTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for projection scheme
 * for 3D MHD physical model.
 *
 * @author Andrea Lani
 */
class MHD3DProjectionEpsVarSet : public MHD3DProjectionVarSet {

public: // classes
  
  typedef MHDProjectionEpsTerm PTERM;
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  MHD3DProjectionEpsVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~MHD3DProjectionEpsVarSet();

  /**
   * Set up the private data and give the maximum size of states
   * physical data to store
   */
  virtual void setup()
  {
    MHD3DProjectionVarSet::setup();
  }
  
  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const = 0;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians() = 0;

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			  RealMatrix& jacobMin,
			  RealVector& eValues,
			  const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"MHD3DProjectionEpsVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"MHD3DProjectionEpsVarSet::computeEigenValuesVectors()");
  }

 
  /**
   * Get the model
   */
  Common::SafePtr<MHDProjectionEpsTerm> getModel() const
  {
    return _modelEps;
  }

 /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, 
				   const RealVector& normal, 
				   RealVector& eValues);
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);
  
protected:

  /// Computes the convective flux projected on a normal
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  /// Computes the physical convective flux
  virtual void computeStateFlux(const RealVector& pdata);
  
private:

  /// acquaintance of the model
  Common::SafePtr<MHDProjectionEpsTerm> _modelEps;
  
}; // end of class MHD3DProjectionEpsVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionEpsVarSet_hh
