#ifndef COOLFluiD_Numerics_FiniteVolume_AtmosphereProps_hh
#define COOLFluiD_Numerics_FiniteVolume_AtmosphereProps_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////

/**
 * This class compute some properties
 *  of the Sun atmosphere
 *
 * @author Alejandro Alvarez
 *
 */
class AtmosphereProps : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  AtmosphereProps(const std::string& name);

  /**
   * Default destructor
   */
  ~AtmosphereProps();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Get Ions Pressure
   */
  CFreal getIonPressure(const CFreal rhoi, const CFreal Ti)
  { return rhoi*Framework::PhysicalConsts::Boltzmann()/1.6726e-27*Ti; }


  /**
   * Get Neutrals Pressure
   */
  CFreal getNeutralPressure(const CFreal rhon, const CFreal Tn)
  { return rhon*Framework::PhysicalConsts::Boltzmann()/1.6726e-27*Tn; }


  /**
  * get the densities and temperatures
  *
  */
  void getDensTemp(const Framework::State* currState, CFreal& rhoi, CFreal& rhon, CFreal& Ti, CFreal& Tn);

  /**
  * Compute chemical reactions frequencies
  *
  */
  void computeChemFreqs(const Framework::State* currState, CFreal& nu_Ion,CFreal& nu_Rec, CFreal& ionsIonizRate, CFreal& neutralsRecombRate);

    /**
  * Compute collisional frequencies
  *
  */
  void computeCollFreqs(const Framework::State* currState, CFreal& nu_in,CFreal& nu_ni, CFreal& nu_en, CFreal& nu_ei);

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
  Framework::DataSocketSource <CFreal> socket_gradPi;

  /// storage of the electric conductivity
  Framework::DataSocketSource <CFreal> socket_gradPn;
  
  /// storage of the x-component of Electric Current
  Framework::DataSocketSource <CFreal> socket_nuIon;
  
  /// storage of the y-component of Electric Current
  Framework::DataSocketSource <CFreal> socket_nuRec;
  
  /// storage of the z-component of Electric Current
  Framework::DataSocketSource <CFreal> socket_nu_in;

  /// storage of Lorenz force
  Framework::DataSocketSource <CFreal> socket_jxB_x;

  /// storage of Lorenz force
  Framework::DataSocketSource <CFreal> socket_jxB_y;

  /// storage of Lorenz force
  Framework::DataSocketSource <CFreal> socket_jxB_z;

  /// storage of Lorenz force
  Framework::DataSocketSource <CFreal> socket_jxB;

  /// storage of the electric current in x
  Framework::DataSocketSource <CFreal> socket_Jx;

  /// storage of the electric current in y
  Framework::DataSocketSource <CFreal> socket_Jy;
  
  /// storage of the electric current in z
  Framework::DataSocketSource <CFreal> socket_Jz;

  /// storage of the magnitude of the electric current
  Framework::DataSocketSource <CFreal> socket_Jtot;
    
  /// pointer to the physical-chemical library In case in the future we take
  /// the cross-sections from library
  //Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceBuilder;
  
  /// gradients of pi
  RealVector m_gradPi;

  /// gradients of pn
  RealVector m_gradPn;

  /// curl of B
  RealVector m_curlB;
  
  /// temporary normal to the face
  RealVector m_normal;  
  
}; // end of class AtmosphereProps

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeMultiFluidMHD_AtmosphereProps_hh



