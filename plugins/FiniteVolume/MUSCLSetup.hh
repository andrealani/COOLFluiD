#ifndef COOLFluiD_Numerics_FiniteVolume_MUSCLSetup_hh
#define COOLFluiD_Numerics_FiniteVolume_MUSCLSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to be executed during the set up
 * of a standard Finite Volume Method.
 *
 * @author Andrea Lani
 */
class MUSCLSetup : public StdSetup {
public:

  /**
   * Constructor.
   */
  explicit MUSCLSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~MUSCLSetup();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Execute Processing actions
   */
  void execute();

private: // data
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSource
  <std::vector<Framework::State*> > socket_stencil;
  
  /// storage for uX
  Framework::DataSocketSource<CFreal> socket_uX;

  /// storage for uY
  Framework::DataSocketSource<CFreal> socket_uY;

  /// storage for uZ
  Framework::DataSocketSource<CFreal> socket_uZ;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MUSCLSetup_hh
