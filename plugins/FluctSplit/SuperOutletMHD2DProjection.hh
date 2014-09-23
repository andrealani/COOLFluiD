#ifndef COOLFluiD_Numerics_FluctSplit_SuperOutletMHD2DProjection_hh
#define COOLFluiD_Numerics_FluctSplit_SuperOutletMHD2DProjection_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic outlet command in 2D
   * for projection scheme
   *
   * @author Radka Keslerova
   *
   */
class SuperOutletMHD2DProjection : public FluctuationSplitCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperOutletMHD2DProjection(const std::string& name);

  /**
   * Default destructor
   */
  ~SuperOutletMHD2DProjection();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on the current TRS
   */
  void executeOnTrs();

protected: // data

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the isUpdated flag
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

private: //data

  /// phi value that is to be fixed
  CFreal _refPhi;

}; // end of class SuperOutlet

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SuperOutletMHD2DProjection_hh
