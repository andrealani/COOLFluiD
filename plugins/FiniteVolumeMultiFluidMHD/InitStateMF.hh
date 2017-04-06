#ifndef COOLFluiD_Numerics_FiniteVolume_InitStateMF_hh
#define COOLFluiD_Numerics_FiniteVolume_InitStateMF_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a initalizing solution command customized for 
   * multi-fluid cases
   *
   * @author Yana Maneva
   *
   */
class InitStateMF : public CellCenterFVMCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit InitStateMF(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~InitStateMF();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

protected:

  /**
   * Execute Processing actions
   */
  virtual void executeOnTrs();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data
  
  /// storage of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
}; // class InitStateMF

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_InitStateMF_hh

