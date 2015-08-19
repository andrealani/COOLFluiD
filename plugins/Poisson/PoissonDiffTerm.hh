#ifndef COOLFluiD_Physics_Poisson_PoissonDiffTerm_hh
#define COOLFluiD_Physics_Poisson_PoissonDiffTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a PoissonDiffTerm.
 *
 * @author Alejandro Alvarez
 *
 */
class PoissonDiffTerm : public Framework::BaseTerm {
  
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
  PoissonDiffTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~PoissonDiffTerm();

  /**
   * Physical data size to be adapted 
   */
  CFuint getDataSize() const
  {
    return 1;
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
    return "PoissonDiffTerm";
  } 
  
private:
// DATA
  
}; // end of class PoissonDiffTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_PoissonDiffTerm_hh
