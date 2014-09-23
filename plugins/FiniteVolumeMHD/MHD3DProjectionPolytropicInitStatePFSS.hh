#ifndef COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionPolytropicInitStatePFSS_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionPolytropicInitStatePFSS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/InitState.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionPolytropicVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents an initializing solution command for MHD3DProjectionPolytropic
   * using PFSS or dipole magnetic field reconstruction and Parker's solar wind with polytropic
   * modelling
   *
   * @author Mehmet Sarp Yalim 
   *
   */
class MHD3DProjectionPolytropicInitStatePFSS : public InitState {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit MHD3DProjectionPolytropicInitStatePFSS(const std::string& name);

  /**
   * Destructor.
   */
  ~MHD3DProjectionPolytropicInitStatePFSS();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Read the PFSS spherical harmonics coefficients from the input file(s)
   */
  void readInputFile();

  /**
   * Compute the Parker's solution for the solar wind velocity and density 
   */
  void computeParkerSolution(const CFreal r, 
			     RealVector& velParkerSpherical,
			     CFreal& rhoParker);
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

protected:

  /**
   * Execute Processing actions
   */
  void executeOnTrs();

  /**
   * Construct the file path from which to read the PFSS spherical harmonics coefficients data
   */
  boost::filesystem::path constructFilename(std::string fileName);

protected: // data

  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD3DProjectionPolytropicVarSet> _varSet;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Alm, to be read from the first file 
  Framework::DataSocketSource<std::vector<CFreal> > socket_AlmrealBegin;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Alm, to be read from the first file
  Framework::DataSocketSource<std::vector<CFreal> > socket_AlmimgBegin;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Blm, to be read from the first file
  Framework::DataSocketSource<std::vector<CFreal> > socket_BlmrealBegin;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Blm, to be read from the first file
  Framework::DataSocketSource<std::vector<CFreal> > socket_BlmimgBegin;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Alm, to be read from the second file 
  Framework::DataSocketSource<std::vector<CFreal> > socket_AlmrealEnd;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Alm, to be read from the second file
  Framework::DataSocketSource<std::vector<CFreal> > socket_AlmimgEnd;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Blm, to be read from the second file
  Framework::DataSocketSource<std::vector<CFreal> > socket_BlmrealEnd;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Blm, to be read from the second file
  Framework::DataSocketSource<std::vector<CFreal> > socket_BlmimgEnd;

  /// socket for the PFSS magnetic field components in Cartesian coordinates to be computed once in the setup phase
  Framework::DataSocketSource<std::vector<CFreal> > socket_BPFSS;

  /// desired accuracy for the Parker solution
  CFreal _epsilon;

  /// radius of the inner boundary
  CFreal _rMin;

  /// radius of the outer boundary
  CFreal _rMax;

  /// Name of the first input file containing PFSS spherical harmonics coefficients 
  std::string _nameBeginPFSSDataFile;

  /// Name of the second input file containing PFSS spherical harmonics coefficients 
  std::string _nameEndPFSSDataFile;

}; // class MHD3DProjectionPolytropicInitStatePFSS

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionPolytropicInitStatePFSS_hh

