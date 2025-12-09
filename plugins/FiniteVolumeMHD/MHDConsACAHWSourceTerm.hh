#ifndef COOLFluiD_Numerics_FiniteVolume_MHDConsACAHWSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHDConsACAHWSourceTerm_hh

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
class MHDConsACAHWSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  MHDConsACAHWSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHDConsACAHWSourceTerm();

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
  Framework::DataSocketSource<CFreal> socket_heating;
  Framework::DataSocketSource<CFreal> socket_radiativeloss;    
  Framework::DataSocketSource<CFreal> socket_wavepressure;
  Framework::DataSocketSource<CFreal> socket_divBCellCenter;
  
//  Framework::DataSocketSource<CFreal> socket_zp;
//  Framework::DataSocketSource<CFreal> socket_zm;
  RealVector _gradP;
  RealVector _gradBx;
  RealVector _gradBy;
  RealVector _gradBz;
  RealVector _gradRho;
  RealVector _gradVx;
  RealVector _gradVy;
  RealVector _gradVz;
  CFint _gravity;
  CFint _PevtsovHeating;
  CFreal _PevtsovHeatingFactor;
  CFint _Manchester;
  CFreal _ManchesterHeatingAmplitude;
  CFreal _ManchesterSigma;
  CFreal _Qh4H_const;
  CFreal _Qlio_AR;
  CFreal _P0_W;
  CFint _alfven_pressure;
  CFint _Qh2_activate;
  CFint _Qh3_activate;
  CFint _Qh4_activate;
  CFint _Qh_lio_activate;
  CFint _divQ;
  CFreal _divQConductivity;
  CFreal _divQalphaCollisionless;
  CFreal _Resistivity;
  CFint _RadiativeLossTerm;
  CFint _wave_pressure;
  CFint _deCompE;
//  CFint _zplus;
//  CFint _zminus;
}; // end of class MHDConsACAHWSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHDConsACAHWSourceTerm_hh
