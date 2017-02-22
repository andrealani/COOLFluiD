#ifndef COOLFluiD_Numerics_FiniteVolume_EnergyIntegration_hh
#define COOLFluiD_Numerics_FiniteVolume_EnergyIntegration_hh

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

  namespace Physics {
    
    namespace MultiFluidMHD {
      class DiffMFMHD2DVarSet;
    }

    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////

/**
 * This class computes the energy variation in the ion+electron cases
 *
 * @author Alejandro Alvarez
 *
 */
class EnergyIntegration : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  EnergyIntegration(const std::string& name);

  /**
   * Default destructor
   */
  ~EnergyIntegration();

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
   * Method writing the solution
   */
  void writeOutputFile();

  /**
   * Method constructing the file
   */
  boost::filesystem::path constructFilename();

  /**
  * Prepare the header of the output file for writing
  */
 void prepareOutputFile(std::ofstream& outputFile);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// save rate
  CFuint  m_saveRate;

  /// Save each iteration to a different Tecplot file (with suffix _iter#).
  bool  m_appendIter;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  bool  m_appendTime;

  /// name of the Output File where to write the electric current
  std::string m_nameOutputFileError;

private: //function

private: //data

  /// physical model var set
  //Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > m_updateVarSet;

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
  
  /// pointer to the physical-chemical library In case in the future we take
  /// the cross-sections from library
  //Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceBuilder;
  
  /// Error total sum over the processors
  CFreal _totalMagneticEner;

  /// Error total sum over the processors
  CFreal _totalElectricEner;

  /// Error total sum over the processors
  CFreal _totalPsiEner;

  /// Error total sum over the processors
  CFreal _totalPhiEner;

  /// Error total sum over the processors
  CFreal _totalKinElecEner;

  /// Error total sum over the processors
  CFreal _totalKinIonsEner;

  /// Error total sum over the processors
  CFreal _totalIntElecEner;

  /// Error total sum over the processors
  CFreal _totalIntIonsEner;

  /// Error total sum over the processors
  CFreal _totalTotalEMEner;
  
  /// Error total sum over the processors
  CFreal _totalTotalElecEner;

  /// Error total sum over the processors
  CFreal _totalTotalIonsEner;

  /// Error total sum over the processors
  CFreal _totalTotalEner;

  /// Error total sum over the processors
  CFuint _totalTotalNbCells;
  
}; // end of class EnergyIntegration

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeMultiFluidMHD_EnergyIntegration_hh



