#ifndef COOLFluiD_Numerics_FiniteVolumeArcJet_CurrentProps_hh
#define COOLFluiD_Numerics_FiniteVolumeArcJet_CurrentProps_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////

/**
 * This class compute some properties of the electric current for the ArcJet module
 *
 * @author Alejandro Alvarez
 *
 */
class CurrentProps : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CurrentProps(const std::string& name);

  /**
   * Default destructor
   */
  ~CurrentProps();

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
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of nstates (states in nodes)
  Framework::DataSocketSink<RealVector> socket_nstates;  
  
  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;  
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;  
  
  /// storage of isOutward
  Framework::DataSocketSink<CFint> socket_isOutward; 
  
  /// storage of normals
  Framework::DataSocketSink<CFreal> socket_normals;  
  
  /// storage of the electric conductivity
  Framework::DataSocketSource <CFreal> socket_sigma;
  
  /// storage of the x-component of Electric Current
  Framework::DataSocketSource <CFreal> socket_Jx;
  
  /// storage of the y-component of Electric Current
  Framework::DataSocketSource <CFreal> socket_Jy;
  
  /// storage of the z-component of Electric Current
  Framework::DataSocketSource <CFreal> socket_Jz;  
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library; 
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceBuilder;
  
  /// gradients of phi
  RealVector m_gradPhi; 
  
  /// temporary normal to the face
  RealVector m_normal;  
  
}; // end of class CurrentProps

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeArcJet_CurrentProps_hh



