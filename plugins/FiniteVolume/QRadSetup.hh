#ifndef COOLFluiD_Numerics_FiniteVolume_QRadSetup_hh
#define COOLFluiD_Numerics_FiniteVolume_QRadSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that creates data needed for radiation 
 * calculations, in particular the array storing the radiative heat source
 *
 * @author Andrea Lani
 */
class QRadSetup : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit QRadSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~QRadSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

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
  
  /// Socket for the states
  Framework::DataSocketSource <CFreal> socket_qrad;
  
  /// Socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
}; // class QRadSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_QRadSetup_hh

