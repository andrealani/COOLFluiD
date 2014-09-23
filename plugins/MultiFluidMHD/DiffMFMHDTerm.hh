#ifndef COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDTerm_hh
#define COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MultiFluidMHDModel.
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez
 *
 */
class DiffMFMHDTerm : public Framework::BaseTerm {
  
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor without arguments
   */
  DiffMFMHDTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~DiffMFMHDTerm();

  /**
   * Physical data size to be adapted 
   */
  CFuint getDataSize() const
  {
    //to test the Braginskii transport
    if (m_nbSpecies == 2) {
      return 10; 				//2 Viscosities + 7 ThermConductiv (ion) + 1 ThermConductiv (neutral)
    }
    else{
      return 2*m_nbSpecies;
    }
  }
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();
   
  /**
   * Get the name
   */
  static std::string getName()
  {
    return "DiffMFMHDTerm";
  }
  
  RealVector& getDynViscosityDim()
  {
    return m_dynViscosityVec;
  }  
  
  RealVector& getThermConductivityDim()
  {
    return m_thermConductivityVec;
  }   
  
  /**
   * @return the number of species
   */
  virtual CFuint getNbSpecies() const {return m_nbSpecies;}   
  
  /**
   * @return the number of species
  */
  virtual bool isBraginskii() const {return m_braginskiiTransport;}
  
  void computeNonInducedEMField(CFreal xCoord, CFreal yCoord);
  
   /**
   * Get the magnetic dipole field and dipole moment values
   */
  RealVector& getNonInducedEMField(CFreal x, CFreal y)
  {
    computeNonInducedEMField(x,y);
    return _NonInducedEMField;
  } 

  
private:
  
  /// dimensional coefficient
  RealVector m_dynViscosityVec;
  
  /// dimensional coefficient
  RealVector m_thermConductivityVec;  
  
  /// number of species
  CFuint m_nbSpecies;
  
  /// Flag to use the Braginskii properties
  bool m_braginskiiTransport;
  
  /// dimensional coefficient to store the options input
  std::vector<CFreal> m_dynViscosity;
  
  /// dimensional coefficient to store the options input
  std::vector<CFreal> m_thermConductivity;  
  
  /// Non Induced part of Electromagnetic field
  RealVector _NonInducedEMField;

  ///introduced non induced electromagnetic field
  std::vector<CFreal> _nonInducedEMField;  
  
}; // end of class DiffMFMHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDTerm_hh
