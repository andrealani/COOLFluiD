#ifndef COOLFluiD_Numerics_FiniteVolume_FilterDiffusion_hh
#define COOLFluiD_Numerics_FiniteVolume_FilterDiffusion_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Environment/FileHandlerOutput.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class CellTrsGeoBuilder;
  }
    
  namespace Numerics {
    
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class flags cells for which diffusive fluxes have to be computed
 *
 * @author Andrea Lani
 *
 */
class FilterDiffusion : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */

  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FilterDiffusion(const std::string& name);

  /**
   * Default destructor
   */
  ~FilterDiffusion();

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
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private: //data
  
  /// flags telling if diffusion terms have to be computed for a given cell
  Framework::DataSocketSink <CFreal> socket_activeDiffusion;
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage for State's
  Framework::DataSocketSink < Framework::State* > socket_gstates;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;  
  
  /// array of physical data
  RealVector m_physicalData;
  
  /// free stream total enthalpy
  CFreal m_freeStreamH;
  
  /// deviation from free stream total enthalpy
  CFreal m_deviationH;
    
}; // end of class FilterDiffusion
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FilterDiffusion_hh
