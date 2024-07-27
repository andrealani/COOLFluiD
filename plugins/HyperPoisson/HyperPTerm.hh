#ifndef COOLFluiD_Physics_HyperPoisson_HyperPTerm_hh
#define COOLFluiD_Physics_HyperPoisson_HyperPTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an HyperP convective
 * physical term
 *
 * @author Rayan Dhib
 * @author Andrea Lani
 */
class HyperPTerm : public Framework::BaseTerm {
public:
    
  #ifdef CF_HAVE_CUDA
  /// nested class defining local options
  template <typename P = NOTYPE>
  class DeviceConfigOptions {
  public:
    /// constructor
    HOST_DEVICE DeviceConfigOptions() {}
    
    /// destructor
    HOST_DEVICE virtual ~DeviceConfigOptions() {}
    
  };
  
  #endif

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between the variable name 
   * and its position in the physical data
   * @pre HyperPTerm::PHI is the computed potential
   */
  enum {PHI=0, BX=1, BY=2, BZ=3};
  
  /**
   * Constructor without arguments
   */
  HyperPTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HyperPTerm();

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
  virtual CFuint getDataSize() const {return 4;}
  
  /**
   * Get the number of scalar vars
   */
  virtual CFuint getNbScalarVars(CFuint i) const {return 0;}
  
  /// Get the number of sets of scalar vars
  virtual CFuint getNbScalarVarSets() const {return 0;}
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Get the name
   */
  static std::string getName() {return "HyperPTerm";}
  
protected:
 

}; // end of class HyperPTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_HyperPoisson_HyperPTerm_hh
