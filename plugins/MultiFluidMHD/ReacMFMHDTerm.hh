#ifndef COOLFluiD_Physics_MultiFluidMHD_ReacMFMHDTerm_hh
#define COOLFluiD_Physics_MultiFluidMHD_ReacMFMHDTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a source term for MultiFluidMHD
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez
 *
 */
class ReacMFMHDTerm : public Framework::BaseTerm {
public:

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   *
   * Tau is the characteristic chemistry time
   */
  enum {TAU=0};

  /**
   * Constructor without arguments
   */
  ReacMFMHDTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ReacMFMHDTerm();

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
    return "ReacMFMHDTerm";
  }
  
}; // end of class ReacMFMHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_ReacMFMHDTerm_hh
