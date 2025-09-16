#ifndef COOLFluiD_FluxReconstructionMethod_PositivityPreservationMHD_hh
#define COOLFluiD_FluxReconstructionMethod_PositivityPreservationMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/BasePhysicality.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MHD {
      class MHD2DProjectionVarSet;
      class MHD3DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that checks and enforces the physicality of 
 * an MHD state, particularly the positivity of the 
 * pressure/density/temperature
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 *
 */
class PositivityPreservationMHD : public BasePhysicality {
public:

  /**
   * Constructor.
   */
  explicit PositivityPreservationMHD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~PositivityPreservationMHD();

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
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: // functions

  /**
   * apply pressure possitivity check
   */
  virtual void enforcePhysicality();
  
  /**
   * Check if the states are physical
   */
  virtual bool checkPhysicality();

protected: // data

  /// minimum allowable value for density
  CFreal m_minDensity;

  /// minimum allowable value for pressure
  CFreal m_minPressure;
  
  /// boolean telling wether to also check the internal solution for physicality
  bool m_checkInternal;
  
  /// boolean to tell whether the complete state is limited or a single variable
  bool m_limCompleteState;
  
  /// physical model 
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> m_varSet;

  /// physical model 2D
  Common::SafePtr<Physics::MHD::MHD2DProjectionVarSet> m_varSet2D;

  /// heat capacity ratio minus one
  CFreal m_gammaMinusOne;
  
  /// variable for physical data of sol
  RealVector m_solPhysData;
  
  /// cell averaged state
  RealVector m_cellAvgState;
  
  /// coefficients for the computation of the cell averaged solution
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;
  
  /// storage of the past states
  //Framework::DataSocketSink< Framework::State*> socket_pastStates;

}; // class PositivityPreservationMHD

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PositivityPreservationMHD_hh
