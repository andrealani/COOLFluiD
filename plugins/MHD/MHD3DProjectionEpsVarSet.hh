#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionEpsVarSet_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionEpsVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
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
class MHD3DProjectionEpsVarSet : public Framework::ConvectiveVarSet {

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
    Framework::ConvectiveVarSet::setup();
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
   * Get the dissipation coefficient for projection scheme
   */
  CFreal getDissipationCoefficient() const
  {
    return getModel()->getDissipationCoefficient();
  }

  /**
   * Get the mass of the external object 
   */
  CFreal getMass() const
  {
    return getModel()->getMass();
  }

  /**
   * Get the reference proton density to non-dimensionalize certain source terms
   */
  CFreal getNRef() const
  {
    return getModel()->getNRef();
  }
  
  /**
   * Get the reference magnetic field to non-dimensionalize certain source terms
   */
  CFreal getBRef() const
  {
    return getModel()->getBRef();
  }
  
  /**
   * Get the reference length to non-dimensionalize certain source terms
   */
  CFreal getLRef() const
  {
    return getModel()->getLRef();
  }
  
  /**
   * Get the reference temperature to non-dimensionalize the equations for the solar wind problem
   */
  CFreal getTRef() const
  {
    return getModel()->getTRef();
  }
  
  /**
   * Get the name of the output file for divB errors
   */
  std::string getNameOutputFile() const
  {
    return getModel()->getNameOutputFile();
  }

  /**
   * Get the frequency of saving the output file for divB errors
   */
  CFuint getOutputFileSaveRate() const
  {
    return getModel()->getOutputFileSaveRate();
  }

  /**
   * Set the transformation matrices between Cartesian and spherical coordinate systems
   */
  void setTransformationMatrices(const RealVector& coords,
		           RealVector& coordsSpherical,
                           RealMatrix& carSphTransMat,
                           RealMatrix& sphCarTransMat);

  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs) 
  {
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3;
  }
  
  /**
   * Get the model
   */
  Common::SafePtr<MHDProjectionEpsTerm> getModel() const
  {
    return _model;
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
  Common::SafePtr<MHDProjectionEpsTerm> _model;
  
}; // end of class MHD3DProjectionEpsVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionEpsVarSet_hh
