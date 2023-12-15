#ifndef COOLFluiD_FluxReconstructionMethod_PhysicalityMFMHD2D_hh
#define COOLFluiD_FluxReconstructionMethod_PhysicalityMFMHD2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/BasePhysicality.hh"

#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MultiFluidMHD {
      //class MultiFluidMHDVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that checks and enforces the physicality of 
 * an Euler/NS 2D state, particularly the positivity of the 
 * pressure/density/temperature
 *
 * @author Firas Ben Ameur
 *
 */
class PhysicalityMFMHD2D : public BasePhysicality {
public:

  /**
   * Constructor.
   */
  explicit PhysicalityMFMHD2D(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~PhysicalityMFMHD2D();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /**
   * apply pressure possitivity check
   */
  virtual void enforcePhysicality();
  
  /**
   * Check if the states are physical: for cons rho and p is checked, for puvt p and T is checked
   */
  virtual bool checkPhysicality();

protected: // data

  /// minimum allowable value for density
  CFreal m_minDensity;

  /// minimum allowable value for pressure
  CFreal m_minPressure;
  
  /// minimum allowable value for temperature
  CFreal m_minTemperature;
  
  /// minimum allowable value for gamma
  CFreal m_minGamma;
  
  /// maximum allowable value for gamma
  CFreal m_maxGamma;
  
  /// Turbulent Variables K, Omega...
  std::vector<CFreal> m_minTurbVars;
  
  /// boolean telling whether to also check the internal solution for physicality
  bool m_checkInternal;
  
  /// boolean to tell whether the complete state is limited or a single variable
  bool m_limCompleteState;
  
  /// boolean telling whether to use the experimental limiter
  bool m_expLim;
  
  /// boolean telling whether to clip gamma
  bool m_clipGamma;
  
  /// physical model 
  //Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet> m_MFMHDVarSet;
  //Common::SafePtr< Framework::MultiScalarTerm< Physics::MultiFluidMHD::EulerMFMHDTerm > > m_MFMHDVarSet;

  /// heat capacity ratio minus one
  CFreal m_gammaMinusOne;
  
  /// variable for physical data of sol
  RealVector m_solPhysData;
  
  /// cell averaged state
  RealVector m_cellAvgState;
  
  /// coefficients for the computation of the cell averaged solution
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;

  /// number of species
  CFuint m_nbSpecies;
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

}; // class PhysicalityMFMHD2D

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PhysicalityMFMHD2D_hh
