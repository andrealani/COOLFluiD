#ifndef COOLFluiD_FluxReconstructionMethod_MHDConsACASourceTerm_hh
#define COOLFluiD_FluxReconstructionMethod_MHDConsACASourceTerm_hh

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
 * This class represents the MHD Source Term for FR
 *
 * @author Ray Vandenhoeck
 *
 */
class MHDConsACASourceTerm : public StdSourceTerm {
  
public:
  
  /**
   * Constructor
   */
  MHDConsACASourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHDConsACASourceTerm();
  
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
  
  virtual void getSToStateJacobian(const CFuint iState);
  
  virtual void getSToGradJacobian(const CFuint iState){};
  
  virtual bool isGradDependent(){return false;};
  
private: //data
  
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal > socket_updateCoeff;
    
  /// socket for gravity values
  Framework::DataSocketSource<CFreal> socket_gravity;
  
  /// storage for Br
  Framework::DataSocketSource<CFreal> socket_Br;
  
  /// storage for Btheta
  Framework::DataSocketSource<CFreal> socket_Btheta;
  
  /// storage for Bphi
  Framework::DataSocketSource<CFreal> socket_Bphi;
  
  /// storage for Vr
  Framework::DataSocketSource<CFreal> socket_Vr;
  
  /// storage for Vtheta
  Framework::DataSocketSource<CFreal> socket_Vtheta;
  
  /// storage for Vphi
  Framework::DataSocketSource<CFreal> socket_Vphi;
    
  bool  m_gravity;
  bool  m_PevtsovHeating;
  CFreal m_PevtsovHeatingFactor;
  bool  m_Manchester;
  CFreal m_ManchesterHeatingAmplitude;
  CFreal m_ManchesterSigma;
  bool  m_divQ;
  CFreal m_divQConductivity;
  CFreal m_divQalphaCollisionless;
  bool  m_ViscosityAndResistivity;
  CFreal m_Viscosity;
  CFreal m_Resistivity;
  bool  m_RadiativeLossTerm;
  bool   m_addUpdateCoeff;
  CFuint m_order;

}; // end of class MHDConsACASourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_MHDConsACASourceTerm_hh
