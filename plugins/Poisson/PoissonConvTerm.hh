#ifndef COOLFluiD_Physics_Poisson_PoissonConvTerm_hh
#define COOLFluiD_Physics_Poisson_PoissonConvTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics { 

    namespace Poisson {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a PoissonConvTerm.
 * 
 * @author Alejandro Alvarez
 *
 */
class PoissonConvTerm : public Framework::BaseTerm {
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
  enum {PHI=0};
  
  /**
   * Constructor without arguments
   */
  PoissonConvTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~PoissonConvTerm();
  
  /**
   * Set physical data
   */
  virtual void setupPhysicalData();
  
  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 1;
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
    return "PoissonConvTerm";
  }
  
protected:
// DATA
  
}; // end of class PoissonConvTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson
 
  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonConvTerm_hh
