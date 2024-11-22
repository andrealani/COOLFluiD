#ifndef COOLFluiD_Physics_PoissonNEQ_PoissonNEQTerm_hh
#define COOLFluiD_Physics_PoissonNEQ_PoissonNEQTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an PoissonNEQ reaction physical term
 *
 * @author Andrea Lani
 */
template <typename BASE>
class PoissonNEQTerm : public BASE {
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
  PoissonNEQTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~PoissonNEQTerm();

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
    return "PoissonNEQTerm";
  }

  /// get flag telling whether to use the axial symmetric model
  bool useAxiModel() const {return m_axiModel;}
  /// Vastalya: get flag telling whether to use the CNEQST(1T) model
  bool useCNEQST() const {return m_CNEQST;}

private:

  /// flag telling whether to use the axial symmetric model
  bool m_axiModel;
  // Vatsalya: flag to see if I want to use NEQST (2T) model or CNEQST (1T). If this is 1, then CNEQST else NEQST
  bool m_CNEQST; 
  
}; // end of class PoissonNEQTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "PoissonNEQ/PoissonNEQTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_PoissonNEQ_PoissonNEQTerm_hh
