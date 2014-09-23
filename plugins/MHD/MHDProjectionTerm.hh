#ifndef COOLFluiD_Physics_MHD_MHDProjectionTerm_hh
#define COOLFluiD_Physics_MHD_MHDProjectionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MHDProjectionTerm.
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
class MHDProjectionTerm : public MHDTerm {
public:
  
  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {PHI=15};

  /**
   * Constructor without arguments
   */
  MHDProjectionTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MHDProjectionTerm();

  /**
   * Configures this object by complementing the
   * implementation in parent.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return MHDTerm::getDataSize() + 1;
  }
  
  /**
   * Get the name
   */
  static std::string getName()
  {
    return "MHDProjectionTerm";
  }

}; // end of class MHDProjectionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionTerm_hh
