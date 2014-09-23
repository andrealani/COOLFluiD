#ifndef COOLFluiD_Physics_KOmega_KOmegaReactionTerm_hh
#define COOLFluiD_Physics_KOmega_KOmegaReactionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an K-Omega source term
 *
 * @author Thomas Wuilbaut
 *
 */
template <typename BASE>
class KOmegaReactionTerm : public BASE {
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
   *
   * The data are:
   * DUMMY - no data for the moment
   */
  // enum {DUMMY=BASE::END};
  
  /**
   * Constructor without arguments
   */
  KOmegaReactionTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~KOmegaReactionTerm();
  
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
    return BASE::getDataSize();
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
    return "KOmegaReactionTerm";
  }

}; // end of class KOmegaReactionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/KOmegaReactionTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_KOmegaReactionTerm_hh
