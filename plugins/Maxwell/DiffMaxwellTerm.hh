#ifndef COOLFluiD_Physics_Maxwell_DiffMaxwellTerm_hh
#define COOLFluiD_Physics_Maxwell_DiffMaxwellTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MaxwellModel.
 *
 * @author Andrea Lani
 *
 */
class DiffMaxwellTerm : public Framework::BaseTerm {

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
  enum {SIGMA=0, END=1};
  
  /**
   * Constructor without arguments
   */
  DiffMaxwellTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffMaxwellTerm();
  
  /**
   * Physical data size
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
    return "DiffMaxwellTerm";
  }
  
}; // end of class DiffMaxwellTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_DiffMaxwellTerm_hh
