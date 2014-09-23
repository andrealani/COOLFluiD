#ifndef COOLFluiD_Physics_ArcJet_ArcJetInductionTerm_hh
#define COOLFluiD_Physics_ArcJet_ArcJetInductionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a ArcJetInductionTerm.
 *
 * @author Andrea Lani
 *
 */
template <typename BASE>
class ArcJetInductionTerm : public BASE {
public:
  
  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {BX=14, BY=15, BZ=16, PHI=17};
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor without arguments
   */
  ArcJetInductionTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ArcJetInductionTerm();

  /**
   * Configures this object by complementing the
   * implementation in parent.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();
 
  /**
   * Get the reference speed for projection scheme
   */
  CFreal getRefSpeed() const {return m_refSpeed;}
  
  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return BASE::getDataSize() + 4;
  }
  
  /**
   * Get the name
   */
  static std::string getName()
  {
    return "ArcJetInductionTerm";
  }
  
private:
  
  /// reference speed
  CFreal m_refSpeed;
  
}; // end of class ArcJetInductionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJet/ArcJetInductionTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetInductionTerm_hh
