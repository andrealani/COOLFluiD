#ifndef COOLFluiD_Numerics_FluctSplit_StrongMirrorAxisymmEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongMirrorAxisymmEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong mirror for axisymmetric Euler2DCons
 *
 * @author Andrea Lani
 *
 */
class StrongMirrorAxisymmEuler2DCons : public FluctuationSplitCom {
public:

  /**
   * Constructor.
   */
  StrongMirrorAxisymmEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongMirrorAxisymmEuler2DCons();

  /**
   * Set up private data
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

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

 }; // end of class StrongMirrorAxisymmEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongMirrorAxisymmEuler2DCons_hh
