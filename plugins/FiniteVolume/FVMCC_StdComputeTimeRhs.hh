#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_StdComputeTimeRhs_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_StdComputeTimeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
/**
 * This class represent a command that computes the time contribution 
 * to the RHS using standard cell center FVM schemes
 *
 * @author Andrea Lani
 *	
 */
class FVMCC_StdComputeTimeRhs : public CellCenterFVMCom {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit FVMCC_StdComputeTimeRhs(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_StdComputeTimeRhs();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> _lss;
    
  /// storage of the States
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
  
  /// storage of the rhs
  Framework::DataSocketSink< CFreal> socket_rhs;
  
  /// array of flags to annull some entries in the matrix diagonal 
  std::vector<bool> _zeroDiagValue;
  
 }; // class FVMCC_StdComputeTimeRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_StdComputeTimeRhs_hh
