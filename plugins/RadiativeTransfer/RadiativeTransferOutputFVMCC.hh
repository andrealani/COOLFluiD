#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferOutputFVMCC_hh
#define COOLFluiD_RadiativeTransfer_RadiativeTransferOutputFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "MathTools/FunctionParser.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Numerics {
    namespace FiniteVolume {
      class CellCenterFVMData;
    }
  }
  
  namespace RadiativeTransfer {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class outputs the radiative wall quantities using Tecplot
 *
 * @author Andrea Lani
 *
 */
class RadiativeTransferOutputFVMCC : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  RadiativeTransferOutputFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~RadiativeTransferOutputFVMCC();

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
  
  /**
   * Compute the values at the wall and write them to file
   */
  virtual void computeWall();
  
  /**
   * Open the Output File and Write the header
   */
  virtual void prepareOutputFileWall();
  
  /**
   * Write the aerodynamic coefficients to file
   */
  virtual void updateOutputFileWall();
  
  /// Update values to be printed
  void updateValuesMat(const CFuint iVar, const CFuint index, const CFreal value)
  {
    cf_assert(iVar < m_valuesMat.nbRows());
    cf_assert(index < m_valuesMat.nbCols());
    m_valuesMat(iVar, index) = value;
  }
  
protected:

  /// builder for Cell Centered FVM TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> m_faceTrsGeoBuilder; 
  
  /// the socket to the data handle of the states
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// the socket to the data handle of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// the socket to the data handle of the nodal states
  Framework::DataSocketSink < RealVector > socket_nstates;
  
  ///  handle to face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// storage of radiative heat flux
  Framework::DataSocketSink<CFreal> socket_qradFluxWall;
  
  /// pointer to the data of the cell centered FVM method
  Common::SafePtr<Numerics::FiniteVolume::CellCenterFVMData> m_fvmccData;
  
  /// mapping between faceIDs and global index
  Common::CFMap<CFuint, CFuint> m_mapTrsFaceToID;
  
  /// local face index
  CFuint m_iFace;
  
  /// current face
  Framework::GeometricEntity* m_currFace;
  
  /// mid face node
  RealVector m_midFaceNode;
  
  /// 2D array storing all values to write to file
  RealMatrix m_valuesMat;
  
  /// list of variable names to write
  std::vector<std::string> m_varNames;
  
  /// name of the surface convergence file
  std::string m_nameOutputFileWall;
  
  ///flag for appending iteration
  bool m_appendIter;
  
  ///flag for appending time
  bool m_appendTime;
  
  /// temperature ID 
  CFuint m_TID;
  
}; /// end of class RadiativeTransferOutputFVMCC

//////////////////////////////////////////////////////////////////////////////

    } /// namespace RadiativeTransfer

} /// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif /// COOLFluiD_RadiativeTransfer_RadiativeTransferOutputFVMCC_hh
