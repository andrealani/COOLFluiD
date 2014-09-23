#ifndef COOLFluiD_Physics_NEQ_NEQReactionTerm_hh
#define COOLFluiD_Physics_NEQ_NEQReactionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a source term for TCNEQ
 *
 * @author Andrea Lani
 *
 */
class NEQReactionTerm : public Framework::BaseTerm {
public:

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   *
   * Tau is the characteristic chemistry time
   */
  enum {TAU=0, END=1};
  
  /**
   * Constructor without arguments
   */
  NEQReactionTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NEQReactionTerm();

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
   * Get the name
   */
  static std::string getName()
  {
    return "NEQReactionTerm";
  }
  
}; // end of class NEQReactionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_NEQReactionTerm_hh
