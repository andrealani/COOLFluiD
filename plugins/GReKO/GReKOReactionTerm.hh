#ifndef COOLFluiD_Physics_GReKO_GReKOReactionTerm_hh
#define COOLFluiD_Physics_GReKO_GReKOReactionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an Gamma-Re_theta-K-Omega source term
 *
 * @author Khalil Bensassi
 *
 */
class GReKOReactionTerm : public Framework::BaseTerm {
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
  enum {DUMMY=0};

  /**
   * Constructor without arguments
   */
  GReKOReactionTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~GReKOReactionTerm();

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
    return "GReKOReactionTerm";
  }

protected:

  /// permeability of the free space
  CFreal m_permeability;

  /// operating frequency of the torch [MHz]
  CFreal m_freqMHz;

}; // end of class GReKOReactionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_GReKOReactionTerm_hh
