#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicVarSet_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "MHD/MHDProjectionPolytropicTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for projection scheme
 * for 3D MHD physical model on the polytropic modelling of the solar wind.
 *
 * @author Mehmet Sarp Yalim
 */
class MHD3DProjectionPolytropicVarSet : public Framework::ConvectiveVarSet {

public: // classes

  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  MHD3DProjectionPolytropicVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~MHD3DProjectionPolytropicVarSet();

  /**
   * Set up the private data and give the maximum size of states
   * physical data to store
   */
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();

    const std::string potentialBType = getModel()->getPotentialBType();

    if (potentialBType == "Dipole") {
       _BDipole.resize(Framework::PhysicalModelStack::getActive()->getDim());
    }

    _tanakaFlux.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
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
    throw Common::NotImplementedException (FromHere(),"MHD3DProjectionPolytropicVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"MHD3DProjectionPolytropicVarSet::computeEigenValuesVectors()");
  }

  /**
   * Compute the convective flux according to Powell, K.G. et. al. (JCP Vol.154 pp.284-309 1999)
   */
  RealVector& computeTanakaFluxPowell99Formulation(const RealVector& pdata,
                                                   const RealVector& normals);

  /**
   * Get the magnetic dipole field and dipole moment values
   */
  RealVector& getMagneticDipole(CFreal x, CFreal y, CFreal z)
  {
    computeMagneticDipole(x,y,z);
    return _BDipole;
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
   * Get the number of l modes to be used in PFSS reconstruction 
   */
  CFuint getNbLModes() const
  {
    return getModel()->getNbLModes();
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
   * Get the x-component of the magnetic dipole moment (mX)
   */
  CFreal getMX() const
  {
    return getModel()->getMX();
  }

  /**
   * Get the y-component of the magnetic dipole moment (mY)
   */
  CFreal getMY() const
  {
    return getModel()->getMY();
  }

  /**
   * Get the z-component of the magnetic dipole moment (mZ)
   */
  CFreal getMZ() const
  {
    return getModel()->getMZ();
  }

  /**
   * Get the potential magnetic field type to model the coronal
   * magnetic field initially: Dipole or PFSS 
   */
  std::string getPotentialBType() const
  {
    return getModel()->getPotentialBType();
  }

  /**
   * Get the name of the method to increase the accuracy in global magnetosphere simulations 
   * especially in the inner magnetosphere: ISLND or Boris 
   */
  std::string getNameAccuracyIncreaserMethod() const
  {
    return getModel()->getNameAccuracyIncreaserMethod();
  }

  /**
   * Compute the magnetic dipole field
   */
  void computeMagneticDipole(CFreal xCoord, CFreal yCoord, CFreal zCoord);

  /**
   * Compute the initial solution for the coronal magnetic field using PFSS technique 
   */
  void computePFSSMagneticField(const RealVector& stateCoordsSpherical,
                                RealVector& BPFSSCartesian,
                                RealMatrix& sphCarTransMat);

  /**
   * Get the reference temperature to non-dimensionalize the equations for the solar wind problem
   */
  CFreal getTRef() const
  {
    return getModel()->getTRef();
  }

  /**
   * Get the ISLND wave speed limit (non-dimensional) 
   */
  CFreal getISLNDLimit() const
  {
    return getModel()->getISLNDLimit();
  }

  /**
   * Get the name of the correction type for projection scheme
   */
  std::string getCorrectionType() const
  {
    return getModel()->getCorrectionType();
  }

  /**
   * Get the name of the output file for divB errors
   */
  std::string getNameOutputFile() const
  {
    return getModel()->getNameOutputFile();
  }

  /**
   * Get the name of the first input file containing the spherical harmonics coefficients for 
   * reconstructing the initial solar coronal magnetic field 
   */
  std::string getNameBeginPFSSCoeffFile() const
  {
    return getModel()->getNameBeginPFSSCoeffFile();
  }

  /**
   * Get the name of the second input file containing the spherical harmonics coefficients for 
   * interpolating the B0 field during the unsteady simulation
   */
  std::string getNameEndPFSSCoeffFile() const
  {
    return getModel()->getNameEndPFSSCoeffFile();
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
  Common::SafePtr<MHDProjectionPolytropicTerm> getModel() const
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
  Common::SafePtr<MHDProjectionPolytropicTerm> _model;

  /// magnetic dipole field
  RealVector _BDipole;

  /// flux vector for Tanaka flux
  RealVector _tanakaFlux;

}; // end of class MHD3DProjectionPolytropicVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicVarSet_hh
