#ifndef COOLFluiD_Physics_MHD_MHDProjectionEpsTerm_hh
#define COOLFluiD_Physics_MHD_MHDProjectionEpsTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHDProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MHDProjectionEpsTerm.
 *
 * @author Andrea Lani
 *
 */
class MHDProjectionEpsTerm : public MHDProjectionTerm {
public:
  
  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {EPSP=16, EPSM=17};
  
  /**
   * Constructor without arguments
   */
  MHDProjectionEpsTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MHDProjectionEpsTerm();
  
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
    return MHDProjectionTerm::getDataSize() + 2;
  }
  
  /**
   * Get the name
   */
  static std::string getName()
  {
    return "MHDProjectionEpsTerm";
  }

}; // end of class MHDProjectionEpsTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionEpsTerm_hh
