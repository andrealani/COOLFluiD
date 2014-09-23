#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_StdComputeTimeRhsCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_StdComputeTimeRhsCoupling_hh

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
class FVMCC_StdComputeTimeRhsCoupling : public CellCenterFVMCom {
public:
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor.
   */
  explicit FVMCC_StdComputeTimeRhsCoupling(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_StdComputeTimeRhsCoupling();
  
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
  
  /// storage of the States
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
  
  /// acquaintance of all the linear system solvers
  std::vector<Common::SafePtr<Framework::LinearSystemSolver> > _lss;
  
  // array of the equation IDs for each LSS
  std::vector<Common::SafePtr<std::vector<CFuint> > > _equations;
  
  // array of jacobian matrices
  std::vector<Common::SafePtr<Framework::LSSMatrix> > _jacobMatrix;
  
  // array of the index mappings from local to global for each LSS
  std::vector<Framework::LSSIdxMapping*> _idxMapping;
  
  /// flag telling if to use global DT (global time stepping)
  bool _useGlobalDT;
  
  /// flag telling if to use analytical transformation matrix
  bool _useAnalyticalMatrix;
  
  /// flag telling to annull the diagonal value in the time contribution 
  /// to the system matrix for a certain subsystem
  std::vector<bool> _annullDiagValue;
  
  /// array of flags to annull some entries in the matrix diagonal 
  std::vector<bool> _zeroDiagValue;
  
 }; // class FVMCC_StdComputeTimeRhsCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_StdComputeTimeRhsCoupling_hh
