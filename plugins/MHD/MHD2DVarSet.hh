#ifndef COOLFluiD_Physics_MHD_MHD2DVarSet_hh
#define COOLFluiD_Physics_MHD_MHD2DVarSet_hh

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
 * This class represents a variable set for a 2D MHD physical model.
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 */
class MHD2DVarSet : public Framework::ConvectiveVarSet {

public: // classes
  
  typedef MHDTerm PTERM;
  
  /**
   * Constructor
   * @see MHDPhysicalModel
   */
  MHD2DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~MHD2DVarSet();

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
    throw Common::NotImplementedException (FromHere(),"MHD2DVarSet::splitJacobian()");
  }

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"MHD2DVarSet::computeEigenValuesVectors()");
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
   * Compute the magnetic dipole field
   */
  void computeMagneticDipole(CFreal xCoord, CFreal yCoord);
  
  /**
   * Compute the convective flux according to Tanaka, T. (JCP Vol.111 pp.381-389 1994)
   */
  RealVector& computeTanakaFlux(const RealVector& data,
				const RealVector& normals);

  /**
   * Compute the convective flux according to Powell, K.G. et. al. (JCP Vol.154 pp.284-309 1999)
   */
  RealVector& computeTanakaFluxPowell99Formulation(const RealVector& data,
		                                   const RealVector& normals)
  {
    throw Common::NotImplementedException (FromHere(),"MHD2DVarSet::computeTanakaFluxPowell99Formulation()");
  }    
  
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

}; // end of class MHD2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DVarSet_hh
