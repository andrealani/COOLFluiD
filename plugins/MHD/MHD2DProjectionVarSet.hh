#ifndef COOLFluiD_Physics_MHD_MHD2DProjectionVarSet_hh
#define COOLFluiD_Physics_MHD_MHD2DProjectionVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "MHDProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for projection scheme
 * for 2D MHD physical model.
 *
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 */
class MHD2DProjectionVarSet : public Framework::ConvectiveVarSet {

public: // classes
  
  typedef MHDProjectionTerm PTERM;
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  MHD2DProjectionVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~MHD2DProjectionVarSet();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();
    _BDipole.resize(Framework::PhysicalModelStack::getActive()->getDim());
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
    throw Common::NotImplementedException (FromHere(),"MHD2DProjectionVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"MHD2DProjectionVarSet::computeEigenValuesVectors()");
  }

  /**
   * Get the dissipation coefficient for projection scheme
   */
  CFreal getDissipationCoefficient() const
  {
    return getModel()->getDissipationCoefficient();
  }

  /**
   * Get the name of the correction type for projection scheme
   */
  std::string getCorrectionType() const
  {
    return getModel()->getCorrectionType();
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
   * Get the ISLND wave speed limit (non-dimensional) 
   */
  CFreal getISLNDLimit() const
  {
    return getModel()->getISLNDLimit();
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
   * Compute the convective flux according to Powell, K.G. et. al. (JCP Vol.154 pp.284-309 1999)
   */
  RealVector& computeTanakaFluxPowell99Formulation(const RealVector& pdata,
                                                   const RealVector& normals);

  /**
   * Compute the magnetic dipole field
   */
  void computeMagneticDipole(CFreal xCoord, CFreal yCoord);
  
  /**
   * Compute the convective flux according to Tanaka, T. (JCP Vol.111 pp.381-389 1994)
   */
  RealVector& computeTanakaFlux(const RealVector& pdata, 
				const RealVector& normals);

  /**
   * Get the magnetic dipole field and dipole moment values
   */
  RealVector& getMagneticDipole(CFreal x, CFreal y)
  {
    computeMagneticDipole(x,y);
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
  Common::SafePtr<MHDProjectionTerm> getModel() const
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
  Common::SafePtr<MHDProjectionTerm> _model;

  /// magnetic dipole field
  RealVector _BDipole;

  /// flux vector for Tanaka flux
  RealVector _tanakaFlux;

}; // end of class MHD2DProjectionVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DProjectionVarSet_hh
