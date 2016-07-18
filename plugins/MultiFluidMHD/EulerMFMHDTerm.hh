#ifndef COOLFluiD_Physics_MultiFluidMHD_EulerMFMHDTerm_hh
#define COOLFluiD_Physics_MultiFluidMHD_EulerMFMHDTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Maxwell/MaxwellProjectionTerm.hh"
#include "Framework/VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a multi-fluid
 * Euler convective
 * physical term
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 *
 */
class EulerMFMHDTerm : public Maxwell::MaxwellProjectionTerm {
  
  enum {START=Maxwell::MaxwellProjectionTerm::END};
  
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {RHO=START,XP=START+1, YP=START+2, ZP=START+3};
  
  /**
   * Constructor without arguments
   */
  EulerMFMHDTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~EulerMFMHDTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Resize the physical data
   */
  virtual void resizePhysicalData(RealVector& physicalData);

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return Maxwell::MaxwellProjectionTerm::getDataSize() + 4;
  }

  /**
   * Get the start of the scalar vars data (mass fractions)
   */
  virtual CFuint getFirstScalarVar(CFuint i) const
  {
    return getDataSize();
  }

  /**
   * Get the number of scalar vars
   */
  virtual CFuint getNbScalarVars(CFuint i) const
  {
    return 0;
  }

  /**
   * Get \f$\gamma\f$
   */
  CFreal getGamma() const
  {
    return _gamma;
  }
   
  /**
   * Get rotational velocity in radiants/sec
   */
  CFreal getOmega() const
  {
    return _omega;
  }
  
  /**
   * Get gas constant in [J/mol K]
   */ 
  CFreal getK() const
  {
    return _K;
  } 
  
  /**
   * Get electronMass in [kg/mol]
   */ 
  CFreal getMolecularMass1() const
  {
    return _molecularMass1;
  } 
  
   /**
   * Get molecularMass2 in [kg/mol]
   */ 
  CFreal getMolecularMass2() const
  {
    return _molecularMass2;
  } 
  
   /**
   * Get molecularMass3 in [kg/mol]
   */ 
  CFreal getMolecularMass3() const
  {
    return _molecularMass3;
  }   
  
  CFreal getPermittivity() const
  {
    return _epsilon;
  } 
  
  CFreal getPermeability() const
  {
    return _mu;
  }
  
  CFreal getLightSpeed() const
  {
    return _lightSpeed;
  }

  bool isLeake() const
  {
    return _isLeake;
  }
  
  /**
   * Get index of pressure
   */
  //static CFuint getPressureTerm()//should be removed?
  //{
  //	return EulerMFMHDTerm::P;
  //}
  
  /**
   * Get dimensional pressure from given pState which is stored in the state
   */
  CFreal getPressureDimFromState(CFreal pState)
  {
    return pState;
  }
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "EulerMFMHDTerm";
  }
  
   /**
   * Compute the Non Induced Electromagnetic field
   */
  void computeNonInducedEMField(CFreal xCoord, CFreal yCoord, CFreal zCoord);

  /**
   * Get the magnetic dipole field and dipole moment values
   */
  RealVector& getNonInducedEMField(CFreal x, CFreal y, CFreal z)
  {
    computeNonInducedEMField(x,y,z);
    return _NonInducedEMField;
  }


protected:

  /// specific heat ratio
  CFreal _gamma;
  
  /// rotational velocity
  CFreal _omega;
   
  /// gas constant [J/ K]
  CFreal _K;

  /// electron mass [kg/mol]
  CFreal _molecularMass1;
  
  /// proton mass [kg/mol]  
  CFreal _molecularMass2;
  
  /// neutral mass(Atomic Hydrogen) [kg/mol]    
  CFreal _molecularMass3;
  
  ///permittivity of free space [F/m]
  CFreal _epsilon;
  
  /// permeabilitty of free space[H/m]
  CFreal _mu;
  
  ///Speed of light
  CFreal _lightSpeed;

  /// checks if an function is used for the electromagnetic field
  bool _useFunction;
  
  /// Non Induced part of Electromagnetic field
  RealVector _NonInducedEMField;
  
  ///introduced non induced electromagnetic field
  std::vector<CFreal> _nonInducedEMField;  

  /// Flag to use Leake's model in the convective term
  bool _isLeake;

  /// Vector for coordinates + time
  RealVector _variables;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

}; // end of class EulerMFMHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_EulerMFMHDTerm_hh
