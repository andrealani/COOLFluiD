#ifndef COOLFluiD_Physics_MHD_MHDProjectionDiffTerm_hh
#define COOLFluiD_Physics_MHD_MHDProjectionDiffTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MHDProjection diffusive model
 *
 * @author Andrea Lani
 *
 */
class MHDProjectionDiffTerm : public Framework::BaseTerm {

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
  enum {MU=0, LAMBDA=1, END=2};
  
  /**
   * Constructor without arguments
   */
  MHDProjectionDiffTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MHDProjectionDiffTerm();
  
  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 2;
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
    return "MHDProjectionDiffTerm";
  }

private:

}; // end of class MHDProjectionDiffTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionDiffTerm_hh
