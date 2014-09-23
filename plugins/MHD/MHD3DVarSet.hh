#ifndef COOLFluiD_Physics_MHD_MHD3DVarSet_hh
#define COOLFluiD_Physics_MHD_MHD3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "MHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 3D MHD physical model.
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 */
class MHD3DVarSet : public Framework::ConvectiveVarSet {

public: // classes
  
  typedef MHDTerm PTERM;
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  MHD3DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~MHD3DVarSet();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();
    _BDipole.resize(PhysicalModelStack::getActive()->getDim());
    _tanakaFlux.resize(PhysicalModelStack::getActive()->getNbEq());
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
    throw Common::NotImplementedException (FromHere(),"MHD3DVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal)
  {
     throw Common::NotImplementedException (FromHere(),"MHD3DVarSet::computeEigenValuesVectors()");
  }

  /**
   * Get the name of the output file for divB errors
   */
  std::string getNameOutputFile() const
  {
    return getModel()->getNameOutputFile();
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
   * Get the frequency of saving the output file for divB errors
   */
  CFuint getOutputFileSaveRate() const
  {
    return getModel()->getOutputFileSaveRate();
  }

  /**
   * Get the mass of the external object 
   */
  CFreal getMass() const
  {
    return getModel()->getMass();
  }

  /**
   * Get the ISLND wave speed limit (non-dimensional) 
   */
  CFreal getISLNDLimit() const
  {
    return getModel()->getISLNDLimit();
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
   * Get the reference temperature to non-dimensionalize the equations for the solar wind problem
   */
  CFreal getTRef() const
  {
    return getModel()->getTRef();
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
   * Compute the magnetic dipole field
   */
  void computeMagneticDipole(CFreal xCoord, CFreal yCoord, CFreal zCoord);

  /**
   * Compute the convective flux according to Tanaka, T. (JCP Vol.111 pp.381-389 1994)
   */
  RealVector& computeTanakaFlux(const RealVector& pdata,
				const RealVector& normals);

  /**
   * Set the transformation matrices between Cartesian and spherical coordinate systems
   */
  void setTransformationMatrices(const RealVector& coords,
                           RealVector& coordsSpherical,
                           RealMatrix& carSphTransMat,
                           RealMatrix& sphCarTransMat);
  
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
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs) 
  {
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3;
  }
  
  /**
   * Get the model
   */
  Common::SafePtr<MHDTerm> getModel() const
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
  Common::SafePtr<MHDTerm> _model;

  /// magnetic dipole field
  RealVector _BDipole;

  /// flux vector for Tanaka flux
  RealVector _tanakaFlux;

}; // end of class MHD3DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DVarSet_hh
