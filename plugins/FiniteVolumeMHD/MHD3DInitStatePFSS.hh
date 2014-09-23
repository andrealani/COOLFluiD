#ifndef COOLFluiD_Numerics_FiniteVolume_MHD3DInitStatePFSS_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD3DInitStatePFSS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/InitState.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents an initializing solution command for MHD3D
   * using PFSS magnetic field reconstruction and Parker's solar wind
   *
   * @author Mehmet Sarp Yalim 
   *
   */
class MHD3DInitStatePFSS : public InitState {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit MHD3DInitStatePFSS(const std::string& name);

  /**
   * Destructor.
   */
  ~MHD3DInitStatePFSS();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Read the PFSS spherical harmonics coefficients from the input file
   */
  void readInputFile();

  /**
   * Compute the initial solution for the coronal magnetic field using PFSS technique 
   */
  void computePFSSMagneticField(const RealVector& stateCoordsSpherical,
				RealVector& BPFSSCartesian,
			        RealMatrix& sphCarTransMat);

  /**
   * Compute the Parker's solution for the solar wind 
   */
  void computeParkerSolution(const CFreal r, 
			     RealVector& velParkerSpherical,
			     CFreal& rhoParker,
			     CFreal& pParker);
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
  Common::SafePtr<Physics::MHD::MHD3DVarSet> _varSet;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Alm, to be used in the boundary conditions
  Framework::DataSocketSource<std::vector<CFreal> > socket_Almreal;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Alm, to be used in the boundary conditions
  Framework::DataSocketSource<std::vector<CFreal> > socket_Almimg;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Blm, to be used in the boundary conditions
  Framework::DataSocketSource<std::vector<CFreal> > socket_Blmreal;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Blm, to be used in the boundary conditions
  Framework::DataSocketSource<std::vector<CFreal> > socket_Blmimg;

  /// mass of the external object to be specified if different than the Sun 
  CFreal _mass;

  /// reference temperature to compute the isothermal speed of sound for Parker's solar wind 
  CFreal _TRef;

  /// desired accuracy for the Parker solution
  CFreal _epsilon;

  /// radius of the inner boundary
  CFreal _rMin;

  /// radius of the outer boundary
  CFreal _rMax;

  /// radius of the source surface for the PFSS model
  CFreal _rSource;

  /// polytropic index
  CFreal _n;

  /// name of the input file containing the PFSS spherical harmonics coefficients data
  std::string _namePFSSDataFile;

}; // class MHD3DInitStatePFSS

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD3DInitStatePFSS_hh

