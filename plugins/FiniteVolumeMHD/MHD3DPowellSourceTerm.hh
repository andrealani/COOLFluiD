#ifndef COOLFluiD_Numerics_FiniteVolume_MHD3DPowellSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD3DPowellSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Physics {
    namespace MHD {
      class MHD3DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents Powell source term for 3D conservative
 * variables
 *
 * @author Mehmet Sarp Yalim
 *
 */
class MHD3DPowellSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  MHD3DPowellSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD3DPowellSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
    
    _globalSockets.createSocketSink<Framework::State*>("states");
  }

  /**
   * Prepare the output file for writing
   */
  void prepareOutputFile(std::ofstream& outputFile);

  /**
   * Write the errors in divB in the output file
   */
  void writeOutputFile();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
 
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
  
protected: // methods

  /**
   * Construct the file path to which to write the RMS source term
   */
  boost::filesystem::path constructFilename();

private: // data

  /// corresponding variable set
  Common::SafePtr<Physics::MHD::MHD3DVarSet> _varSet;

  /// socket for divB values at the cell centers to pass to FVMCC_ComputeRHS for weighted averaging
  Framework::DataSocketSource<CFreal> socket_divBCellCenter;

  /// socket for average B values in x-direction at the cell faces to be used in the computation of divB in Powell source term
  Framework::DataSocketSink<CFreal> socket_avgBxFace;
       
  /// socket for average B values in y-direction at the cell faces to be used in the computation of divB in Powell source term
  Framework::DataSocketSink<CFreal> socket_avgByFace;
           
  /// socket for average B values in z-direction at the cell faces to be used in the computation of divB in Powell source term
  Framework::DataSocketSink<CFreal> socket_avgBzFace;
  
  /// magnetic dipole field vector
  RealVector _BDipole;

  /// Array storing the divB error values in the domain
  RealVector _divB;

  /// Maximum divB error in the whole domain
  CFreal _divBMax;

  /// Minimum divB error in the whole domain
  CFreal _divBMin;

  /// Storage for choosing when to save the divB error file
  CFuint _saveRate;

  /// Name of output file where to write divB errors
  std::string _nameOutputFile;

  /// MHD physical data
  RealVector _physicalData;

  /// MHD physical data
  RealVector _dataLeftState;

   /// MHD physical data
  RealVector _dataRightState;

}; // end of class MHD3DPowellSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD3DPowellSourceTerm_hh
