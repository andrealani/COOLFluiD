#ifndef COOLFluiD_Numerics_FiniteVolume_MHDConsACASourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHDConsACASourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents source term for artificial compressibility
 * method 
 *
 * @author Mehmet Sarp Yalim
 *
 */
class MHDConsACASourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
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
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

 /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private:
  
  /// socket for gravity values
  Framework::DataSocketSource<CFreal> socket_gravity;
    
  CFint _gravity;
  CFint _PevtsovHeating;
  CFreal _PevtsovHeatingFactor;
  CFint _Manchester;
  CFreal _ManchesterHeatingAmplitude;
  CFreal _ManchesterSigma;
  CFint _divQ;
  CFreal _divQConductivity;
  CFreal _divQalphaCollisionless;
  CFint _ViscosityAndResistivity;
  CFreal _Viscosity;
  CFreal _Resistivity;
  CFint _RadiativeLossTerm;

}; // end of class MHDConsACASourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHDConsACASourceTerm_hh
