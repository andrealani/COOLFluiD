#ifndef COOLFluiD_Physics_Maxwell_MaxwellProjectionAdimTerm_hh
#define COOLFluiD_Physics_Maxwell_MaxwellProjectionAdimTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConvMaxwellTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MaxwellProjectionAdimTerm.
 *
 * @author Alejandro Alvarez Laguna
 * @author Andrea Lani
 *
 */
class MaxwellProjectionAdimTerm : public ConvMaxwellTerm {
  
  enum {START=ConvMaxwellTerm::END}; 

public:
    
  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data. PSI for divB correction and 
   * PHI for divE correction
   */
  enum {PSI = START, PHI = START + 1, END = START + 2};

  /**
   * Constructor without arguments
   */
  MaxwellProjectionAdimTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MaxwellProjectionAdimTerm();

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
    return ConvMaxwellTerm::getDataSize() + 2;
  }

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "MaxwellProjectionAdimTerm";
  }
  
}; // end of class MaxwellProjectionAdimTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_MaxwellProjectionAdimTerm_hh
