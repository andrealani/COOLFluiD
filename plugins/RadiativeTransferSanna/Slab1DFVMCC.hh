#ifndef COOLFluiD_RadiativeTransfer_Slab1DFVMCC_hh
#define COOLFluiD_RadiativeTransfer_Slab1DFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "MathTools/FunctionParser.hh"
#include "MathTools/RealVector.hh" 
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class RadiationLibrary;
  }
  
  namespace Numerics {
    namespace FiniteVolume {
      class CellCenterFVMData;
    }
  }
  
  namespace RadiativeTransfer {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Wall values and aerodynamic coefficients
 *
 * @author Andrea Lani
 * @author Alessandro Munafo'
 *
 */
class Slab1DFVMCC : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Slab1DFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~Slab1DFVMCC();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();

protected:

  /**
   * Execute on a set of dofs
   */
  virtual void executeOnTrs();
  
  /// compute only on the stagnation line
  bool computeOnStagnationLine() const
  {
    return (m_stagnationPointXYZ.size() > 0);
  } 
  
  /// compute the radiative heat flux on the wall
  void computeWallRadHeatFlux();
  
protected:
 
  /// physical properties library
  Common::SelfRegistPtr<Framework::RadiationLibrary> m_radLibrary;
  
  /// the socket to the radiative heat as coming out from the line integration 
  Framework::DataSocketSource < CFreal > socket_qrad;
  
  /// the socket to the radiative heat flux at the wall faces
  Framework::DataSocketSource < CFreal > socket_qradFluxWall;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// stagnation point coordinate vector
  RealVector m_stagPoint;
  
  /// total number of wall faces
  CFuint m_nbWallFaces;
  
  /// mapping between nodeID and stateID
  std::vector<CFuint> m_nodeIdToStateId;
  
  /// local (to the TRS) IDs of the TRS faces ordered according to x direction
  std::vector<CFuint> m_orderedFaceIDs;
  
  /// face ID on the stagnation point
  Common::SafePtr<std::vector<CFuint> > m_stagnationLineCells;
    
  /// pointer to the data of the cell centered FVM method
  Common::SafePtr<Numerics::FiniteVolume::CellCenterFVMData> m_fvmccData;
  
  /// array of cell IDs obtained while marching from wall to outer boundary
  std::vector<std::vector<CFuint>* > m_meshByLines;
  
  /// distances between face centers along the line  
  std::vector<std::vector<CFreal>* > m_deltal;
  
  /// inner TRS boundary names
  std::string m_radLibraryName;
  
  /// user defined stagnation point coordinates
  std::vector<CFreal> m_stagnationPointXYZ;
  
}; /// end of class Slab1DFVMCC

//////////////////////////////////////////////////////////////////////////////

    } /// namespace RadiativeTransfer

} /// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif /// COOLFluiD_RadiativeTransfer_Slab1DFVMCC_hh
