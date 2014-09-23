#ifndef COOLFluiD_Physics_NavierStokes_NSTurbTerm_hh
#define COOLFluiD_Physics_NavierStokes_NSTurbTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "NSTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokes with Turbulence Model.
 *
 * @author Thomas Wuilbaut
 *
 */
class NSTurbTerm : public NSTerm {
public:

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {MUT=3};

  /**
   * Constructor without arguments
   */
  NSTurbTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NSTurbTerm();

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return  NSTerm::getDataSize() + 1;
  }

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "NSTurbTerm";
  }


}; // end of class NSTurbTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NSTurbTerm_hh
