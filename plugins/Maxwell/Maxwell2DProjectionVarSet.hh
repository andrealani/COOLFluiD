#ifndef COOLFluiD_Physics_Maxwell_Maxwell2DProjectionVarSet_hh
#define COOLFluiD_Physics_Maxwell_Maxwell2DProjectionVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "MaxwellVarSet.hh"
#include "Maxwell/MaxwellProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for projection scheme
 * for 2D Maxwell physical model.
 *
 * @author Alejandro Alvarez Laguna
 * @author Andrea Lani
 */
class Maxwell2DProjectionVarSet : public MaxwellVarSet {
public: // classes
  
  typedef MaxwellProjectionTerm PTERM;
  typedef Maxwell2DProjectionVarSet EULERSET;
  
  /**
   * Constructor
   * @see MaxwellModel
   */
  Maxwell2DProjectionVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Maxwell2DProjectionVarSet();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();
    
    // set EquationSetData
    Maxwell2DProjectionVarSet::getEqSetData().resize(1);
    Maxwell2DProjectionVarSet::getEqSetData()[0].setup(0,0,8);   
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
			  const RealVector& normal) = 0;

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				     RealMatrix& leftEv,
				     RealVector& eValues,
				     const RealVector& normal) = 0;
  
  
  /**
   * Get the model
   */
  Common::SafePtr<MaxwellProjectionTerm> getModel() const
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

   /// @returns the number of equations of this VarSet
  CFuint getNbEqs() const { return 8; }
  
private:

  /// acquaintance of the model
  Common::SafePtr<MaxwellProjectionTerm> _model;

}; // end of class Maxwell2DProjectionVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_Maxwell2DProjectionVarSet_hh
