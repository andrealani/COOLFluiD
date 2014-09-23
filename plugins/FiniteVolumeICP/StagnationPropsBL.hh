#ifndef COOLFluiD_Numerics_Finite_VolumeICP_StagnationPropsBLComm_hh
#define COOLFluiD_Numerics_Finite_VolumeICP_StagnationPropsBLComm_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////

/**
 * This class implements the Lorentz force for the induction equations
 *
 * @author Vladimir Prokop
 * @author Emanuele Sartori
 *
 */
class StagnationPropsBL : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  StagnationPropsBL(const std::string& name);

  /**
   * Default destructor
   */
  ~StagnationPropsBL();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

/**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: //function

private: //data

  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// location of the max X coordinate found and used to transform
  CFreal m_xmax;

  /// wall closest states, ordered and x is transformed such that zero is the nose of the probe
  std::vector<std::pair<CFreal,CFuint> > m_localstagnationline;

  /// wall closest states, ordered and x is transformed such that zero is the nose of the probe, only valid on process zero
  /// this is the globally collected array
  std::vector<std::pair<CFreal,CFuint> > m_globalstagnationline;

  /// receive sizes for gathering for process zero, only valid on process zero
  std::vector<int> m_nrcv;

  /// receive displacements for gathering for process zero, only valid on process zero
  std::vector<int> m_dspl;

  /// X coordinate where the torch ends and the chamber starts
  CFreal m_torchexitxcoord;

}; // end of class StagnationPropsBL

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeICP_StagnationPropsBLComm_hh



