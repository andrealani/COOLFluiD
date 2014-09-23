#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHS_LES_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHS_LES_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FiniteVolume/FVMCC_ComputeRHS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace LES {
    class LESVarSet;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes but modified for use with LES
 *
 * @author Willem Deconinck
 *
 */
class FVMCC_ComputeRHS_LES : public FVMCC_ComputeRHS {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRHS_LES(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRHS_LES();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Un Setup private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

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
  
  Common::SafePtr<LES::LESVarSet> m_lesVarSet;
  
}; // class FVMCC_ComputeRHS_LES

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHS_LES_hh

