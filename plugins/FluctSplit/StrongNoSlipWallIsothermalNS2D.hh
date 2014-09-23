#ifndef COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalNS2D_hh
#define COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalNS2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for Euler2D
 *
 * @author Andrea Lani
 */
class StrongNoSlipWallIsothermalNS2D : public FluctuationSplitCom {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  StrongNoSlipWallIsothermalNS2D(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongNoSlipWallIsothermalNS2D();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Set up private data
   */
  void setup();
  
protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// flag telling if initialization is performed
  bool _useForInitialization;
  
  /// dimensional wall temperature
  CFreal _TWall;
  
}; // end of class StrongNoSlipWallIsothermalNS2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalNS2D_hh
