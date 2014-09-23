#ifndef COOLFluiD_Physics_Maxwell_MaxwellAdimTerm_hh
#define COOLFluiD_Physics_Maxwell_MaxwellAdimTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConvMaxwellTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MaxwellAdimTerm.
 *
 * @author Alejandro Alvarez Laguna
 * @author Andrea Lani
 *
 */
class MaxwellAdimTerm : public ConvMaxwellTerm {
  
public:

  /**
   * Constructor without arguments
   */
  MaxwellAdimTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MaxwellAdimTerm();

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
    return ConvMaxwellTerm::getDataSize();
  }

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "MaxwellAdimTerm";
  }
  
}; // end of class MaxwellAdimTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_MaxwellAdimTerm_hh
