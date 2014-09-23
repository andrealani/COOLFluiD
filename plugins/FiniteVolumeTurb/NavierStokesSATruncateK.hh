#ifndef COOLFluiD_Numerics_FiniteVolume_NavierStokesSATruncateK_hh
#define COOLFluiD_Numerics_FiniteVolume_NavierStokesSATruncateK_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class is a post processing command that truncates the entry RhoK in the states.
 *
 * @author Tiago Quintino
 * @author Joao Pinto
 *
 */
class NavierStokesSATruncateK : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesSATruncateK(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesSATruncateK();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

private: //data

  /// socket for Node's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // end of class NavierStokesSATruncateK

//////////////////////////////////////////////////////////////////////////////

  } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokesSATruncateK_hh
