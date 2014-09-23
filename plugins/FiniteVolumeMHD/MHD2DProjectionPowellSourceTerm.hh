#ifndef COOLFluiD_Numerics_FiniteVolume_MHD2DProjectionPowellSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD2DProjectionPowellSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Physics {
    namespace MHD {
      class MHD2DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents Powell source term for projection scheme
 * for 2D conservative variables
 *
 * @author Mehmet Sarp Yalim
 *
 */
class MHD2DProjectionPowellSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  MHD2DProjectionPowellSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD2DProjectionPowellSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);

    _sockets.createSocketSink<Framework::State*, Framework::GLOBAL>("states");
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
   * Set the variable set
   * @pre the input pointer is non const to allow dynamic_cast
   */
  void setVarSet(Common::SafePtr<Framework::ConvectiveVarSet> varSet);

  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
         RealVector& source);

protected: // methods

  /**
   * Construct the file path to which to write the RMS source term
   */
  boost::filesystem::path constructFilename();

private: // data

  /// corresponding variable set
  Common::SafePtr<Physics::MHD::MHD2DProjectionVarSet> _varSet;

  /// Array storing the divB error values in the domain
  RealVector _divB;

  /// Maximum divB error in the whole domain
  CFreal _divBMax;

  /// Minimum divB error in the whole domain
  CFreal _divBMin;

  /// left eigen values
  RealVector _leftEvals;

  /// right eigen values
  RealVector _rightEvals;

  /// unit normal
  RealVector _unitNormal;

  /// MHD physical data
  RealVector _physicalData;

  /// MHD physical data
  RealVector _dataLeftState;

   /// MHD physical data
  RealVector _dataRightState;

  /// Storage for choosing when to save the divB error file
  CFuint _saveRate;

  /// Name of output file where to write divB errors
  std::string _nameOutputFile;

}; // end of class MHD2DProjectionPowellSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD2DProjectionPowellSourceTerm_hh
