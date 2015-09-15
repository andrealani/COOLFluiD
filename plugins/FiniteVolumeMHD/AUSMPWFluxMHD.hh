#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPWFluxMHD_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPWFluxMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeMHD/AUSMFluxMHD.hh"
#include "Framework/EquationSetData.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MeshData.hh"

#include "MHD/MHDTerm.hh"
#include "Maxwell/Maxwell.hh"

#include "MHD/MHDProjectionTerm.hh"
#include "Maxwell/MaxwellProjection.hh"
#include "Maxwell/MaxwellProjectionAdim.hh"
#include "Maxwell/MaxwellProjectionTerm.hh"
#include "Maxwell/MaxwellProjectionAdimTerm.hh"

#include "MHD/MHD3DProjectionVarSet.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"

#include "MHD/MHD3DVarSet.hh"
#include "MHD/MHD2DVarSet.hh"
#include "Maxwell/Maxwell3DVarSet.hh"
#include "Maxwell/Maxwell2DVarSet.hh"
#include "Maxwell/MaxwellVarSet.hh"

#include "MHD/MHD2DCons.hh"
#include "MHD/MHD3DCons.hh"
#include "Maxwell/Maxwell2DCons.hh"
#include "Maxwell/Maxwell3DCons.hh"

#include "FiniteVolume/RoeFlux.hh"
#include "FiniteVolume/LaxFriedFlux.hh"
#include "Environment/ModuleRegister.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM flux
 *
 * @author Alejandro Alvarez Laguna
 * @author Jean-Cedric Chkair
 */
template <class UPDATEVAR>
class AUSMPWFluxMHD : public AUSMFluxMHD<UPDATEVAR> {

public:
  
  /**
   * Constructor
   */
  AUSMPWFluxMHD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPWFluxMHD();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

//  /**
//   * Set up private data
//   */
//  virtual void setup();


protected:

  /**
   * Compute the interface Mach
   */
  virtual void computeInterfaceMach();

//  /**
//   * Compute the interface pressure flux
//   */
//  virtual void computePressureFlux();

  /**
   * Compute the functions of references
   */
  virtual void computeReferences();

  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFluxImpl(RealVector& result);

private:
  
  CFreal m_beta;					// beta coeficient
  CFreal m_alpha;					// alpha coeficient

protected:
  
  /// mach number 4Plus
  CFreal m_M4Plus;

  /// mach number 4Minus
  CFreal m_M4Minus;
  
  /// mach number M5Plus
  CFreal m_M5Plus;
  
  /// mach number M5Minus
  CFreal m_M5Minus;

  /// pressure P5Plus
  CFreal m_P5Plus;

  /// pressure P5Minus
  CFreal m_P5Minus;

  /// normal component of left magnetic field
  CFreal m_BnL;

  /// normal component of right magnetic field
  CFreal m_BnR;

  /// average of normal magnetic fields
  CFreal m_Bn;

  /// left magnetic field square
  CFreal B2L;

  /// right magnetic field square
  CFreal B2R;

  /// scalar velocity-magnetic field of left side
  CFreal VdotBL;

  /// scalar velocity-magnetic field of right side
  CFreal VdotBR;

  /// left velocity square
  CFreal V2L;

  /// right velocity square
  CFreal V2R;

  /// left magnetic field square
  CFreal m_EL;

  /// right magnetic field square
  CFreal m_ER;

  /// left total pressure
  CFreal m_pgL;

  /// right total pressure
  CFreal m_pgR;

}; // end of class AUSMPWFluxMHD

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPWFluxMHD.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPWFluxMHD_hh
