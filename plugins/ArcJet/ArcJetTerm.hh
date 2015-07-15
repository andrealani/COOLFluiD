#ifndef COOLFluiD_Physics_ArcJet_ArcJetTerm_hh
#define COOLFluiD_Physics_ArcJet_ArcJetTerm_hh

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
class ArcJetTerm : public BASE {
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
   */
  enum {SIGMA=BASE::END};
  
  /**
   * Constructor without arguments
   */
  ArcJetTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ArcJetTerm();

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
    return BASE::getDataSize() + 1;
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
    return "ArcJetTerm";
  }
    
}; // end of class ArcJetTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJet/ArcJetTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetTerm_hh
