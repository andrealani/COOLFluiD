#ifndef COOLFluiD_Physics_Maxwell_Maxwell3DVarSetAdim_hh
#define COOLFluiD_Physics_Maxwell_Maxwell3DVarSetAdim_hh

//////////////////////////////////////////////////////////////////////////////

#include "MaxwellVarSet.hh"
#include "Maxwell/MaxwellAdimTerm.hh"
#include "Framework/EquationSetData.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 3D Maxwell physical model for Adimensional equations.
 *
 * @author Alejandro Alvarez Laguna
 * @author Andrea Lani
 */
class Maxwell3DVarSetAdim : public MaxwellVarSet {
public: // classes
  
  typedef Maxwell3DVarSetAdim EULERSET;
  
  /**
   * Constructor
   * @see MaxwellPhysicalModel
   */
  Maxwell3DVarSetAdim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Maxwell3DVarSetAdim();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup()
  {
    std::cout << "Maxwell3DVarSetAdim::setup()/n";
    Framework::ConvectiveVarSet::setup();
    
    // set EquationSetData
    Maxwell3DVarSetAdim::getEqSetData().resize(1);
    Maxwell3DVarSetAdim::getEqSetData()[0].setup(0,0,6);   
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
   * Get some data corresponding to the subset of equations related with
   * this variable set
   * @pre The most concrete ConvectiveVarSet will have to set these data
   */
  static std::vector<Framework::EquationSetData>& getEqSetData()
  {
    static std::vector<Framework::EquationSetData> eqSetData;
    return eqSetData;
  }

  /**
   * Get the model
   */
  Common::SafePtr<MaxwellAdimTerm> getModel() const
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
  CFuint getNbEqs() const { return 6; std::cout << "getNbEqs()\n"; }
  
private:

  /// acquaintance of the model
  Common::SafePtr<MaxwellAdimTerm> _model;

}; // end of class Maxwell3DVarSetAdim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_Maxwell3DVarSetAdim_hh
