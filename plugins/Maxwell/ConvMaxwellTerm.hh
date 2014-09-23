#ifndef COOLFluiD_Physics_Maxwell_ConvMaxwellTerm_hh
#define COOLFluiD_Physics_Maxwell_ConvMaxwellTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an Euler convective
 * physical term
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 *
 */
class ConvMaxwellTerm : public Framework::BaseTerm {
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
  enum {BX=0, BY=1, BZ=2, EX=3, EY=4, EZ=5, END=6};
  //enum {BX=0, BY=1, BZ=2, EX=3, EY=4, EZ=5, PSI=6, END=7};
  
  /**
   * Constructor without arguments
   */
  ConvMaxwellTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConvMaxwellTerm();

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
    return 6;
    //return 7;
  }
        
//   /**
//    * Get the name of the correction type for projection scheme
//    */
//   std::string getCorrectionType() const
//   {
//     return _correctionType;
//   }  

  /**
   * Get the start of the scalar vars data (mass fractions)
   */
  virtual CFuint getFirstScalarVar(CFuint i) const
  {
    return 6;
    //return 7;
  }

  /**
   * Get the number of scalar vars
   */
  virtual CFuint getNbScalarVars(CFuint i) const
  {
    return 0;
  }
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );  

  /**
   * Get DivB Cleaning Technique Constant 
   */  
  
  CFreal getDivBCleaningConst() const 
  {
    return _divBCleaningConst;
  }  
  
  /**
   * Get DivE Cleaning Technique Constant 
   */  
  
  CFreal getDivECleaningConst() const 
  {
    return _divECleaningConst;
  }   
  
  CFreal getDivBAdimCleaningConst() const 
  {
    return _divBAdimCleaningConst;
  }      
  
  /**
   * Get LightSpeed 
   */
  
  CFreal getLightSpeed() const {return _LightSpeed;}

  /**
   * Get electronMass 
   */
  
  CFreal getElectronMass() const {return _electronMass;}  
   
  /**
   * Get electronCharge 
   */
  
  CFreal getElectronCharge() const {return _electronCharge;}    
  
  /**
   * Get protonMass 
   */
  
  CFreal getProtonMass() const {return _protonMass;}     
  
   /**
   * Get neutralMass 
   */
  
  CFreal getNeutralMass() const {return _neutralMass;} 
  
  /**
   * Get solarGravity 
   */
  
  CFreal getSolarGravity() const {return _solarGravity;}  
  
  /**
   * Get the name
   */
  
  static std::string getName()
  {
    return "ConvMaxwellTerm";
  }
  
//   /**
//    * Get the name of the output file for divB errors
//    */
//   std::string getNameOutputFile() const 
//   {
//     return _nameOutputFile;
//   }

  
//   /**
//    * Get the frequency of saving the output file for divB errors
//    */
//   CFuint getOutputFileSaveRate() const 
//   {
//     return _saveRate;
//   }

protected:
  
//   /// Storage for choosing when to save the divB error output file
//   CFuint _saveRate;

//   /// Name of Output File where to write the coefficients
//   std::string _nameOutputFile;



//  /// Name of correction for projection scheme
//  std::string 	    _correctionType;
  
  /// adimensional parameter necessary for the hyperbolic 
  ///divB cleaning technique
  
  CFreal            _divBCleaningConst;
  
  /// adimensional parameter necessary for the hyperbolic 
  ///divE cleaning technique
  CFreal            _divECleaningConst;  
  
  /// adimensional parameter necessary for the hyperbolic 
  ///divB cleaning technique
  CFreal            _divBAdimCleaningConst;   
  
  CFreal            _LightSpeed;
  CFreal            _electronMass;
  CFreal	    _electronCharge;
  CFreal	    _protonMass;
  CFreal	    _neutralMass;
  CFreal	    _solarGravity;

}; // end of class ConvMaxwellTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_ConvMaxwellTerm_hh
