#ifndef COOLFluiD_Physics_GammaAlpha_GammaAlphaReactionTerm_hh
#define COOLFluiD_Physics_GammaAlpha_GammaAlphaReactionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a gamma-alpha source term
 *
 * @author Ray Vandenhoeck
 *
 */
class GammaAlphaReactionTerm : public Framework::BaseTerm {
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
  GammaAlphaReactionTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~GammaAlphaReactionTerm();

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
    return "GammaAlphaReactionTerm";
  }

protected:

  /// permeability of the free space
  CFreal m_permeability;

  /// operating frequency of the torch [MHz]
  CFreal m_freqMHz;

}; // end of class GammaAlphaReactionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_GammaAlphaReactionTerm_hh
