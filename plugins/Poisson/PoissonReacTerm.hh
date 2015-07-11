#ifndef COOLFluiD_Physics_Poisson_PoissonReacTerm_hh
#define COOLFluiD_Physics_Poisson_PoissonReacTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an Poisson reaction physical term
 *
 * @author Alejandro Alvarez
 */

class PoissonReacTerm : public Framework::BaseTerm {
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
   * PHI - electric conductivity
   * SIGMA - permeability of free space
   */
  enum {PHI=0, SIGMA=1};

  /**
   * Constructor without arguments
   */
  PoissonReacTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~PoissonReacTerm();

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
    return 2;
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
    return "PoissonReacTerm";
  }

protected:
  //DATA

}; // end of class PoissonReacTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonReacTerm_hh
