#ifndef COOLFluiD_Physics_ArcJet_ArcJetReactionTerm_hh
#define COOLFluiD_Physics_ArcJet_ArcJetReactionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an ArcJet reaction physical term
 *
 * @author Andrea Lani
 */
template <typename BASE>
class ArcJetReactionTerm : public BASE {
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
   * SIGMA - electric conductivity
   * MU0 - permeability of free space
   */
  enum {SIGMA=0, MU0=1};

  /**
   * Constructor without arguments
   */
  ArcJetReactionTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ArcJetReactionTerm();

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
    return BASE::getDataSize() + 2;
  }
  
  /**
   * Get the permeability of free space
   */
  CFreal getPermeability() const
  {
    return m_permeability;
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
    return "ArcJetReactionTerm";
  }
  
protected:
  
  /// permeability of the free space
  CFreal m_permeability;

}; // end of class ArcJetReactionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJet/ArcJetReactionTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetReactionTerm_hh
