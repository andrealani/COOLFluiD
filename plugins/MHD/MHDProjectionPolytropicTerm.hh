#ifndef COOLFluiD_Physics_MHD_MHDProjectionPolytropicTerm_hh
#define COOLFluiD_Physics_MHD_MHDProjectionPolytropicTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHDProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MHDProjectionPolytropicTerm.
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
class MHDProjectionPolytropicTerm : public MHDProjectionTerm {

public:

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {N=16};

  /**
   * Constructor without arguments
   */
  MHDProjectionPolytropicTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MHDProjectionPolytropicTerm();

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
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return MHDProjectionTerm::getDataSize() + 1;
  }

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "MHDProjectionPolytropicTerm";
  }

}; // end of class MHDProjectionPolytropicTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionPolytropicTerm_hh
