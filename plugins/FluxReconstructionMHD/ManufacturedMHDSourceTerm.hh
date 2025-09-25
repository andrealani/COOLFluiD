#ifndef COOLFluiD_FluxReconstructionMethod_ManufacturedMHDSourceTerm_hh
#define COOLFluiD_FluxReconstructionMethod_ManufacturedMHDSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Common/SafePtr.hh"
#include "FluxReconstructionMethod/StdSourceTerm.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Source Term for the MHD manufactured case for FR
 *
 * @author Rayan Dhib
 *
 */
class ManufacturedMHDSourceTerm : public StdSourceTerm {
  
public:
  
  /**
   * Constructor
   */
  ManufacturedMHDSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ManufacturedMHDSourceTerm();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * add the source term
   */
  void addSourceTerm(RealVector& resUpdates);
  
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
  
  //virtual void getSToStateJacobian(const CFuint iState);
  
  virtual void getSToGradJacobian(const CFuint iState){};
  
  virtual bool isGradDependent(){return false;};
  
protected:

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );


private: //data
  
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal > socket_updateCoeff;
    
  bool   m_addUpdateCoeff;
  CFuint m_order;

}; // end of class ManufacturedMHDSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ManufacturedMHDSourceTerm_hh
