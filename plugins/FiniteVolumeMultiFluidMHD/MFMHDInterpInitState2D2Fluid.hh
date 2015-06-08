#ifndef COOLFluiD_Numerics_FiniteVolume_MFMHDInterpInitState2D2Fluid_hh
#define COOLFluiD_Numerics_FiniteVolume_MFMHDInterpInitState2D2Fluid_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/InitState.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * Interpolation from a table to the grid of initial state for 2Fluid in 2D grids
   *
   * @author Alejandro Alvarez
   * @author Yana Maneva
   *
   */
class MFMHDInterpInitState2D2Fluid : public InitState {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit MFMHDInterpInitState2D2Fluid(const std::string& name);

  /**
   * Destructor.
   */
  ~MFMHDInterpInitState2D2Fluid();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Read the Input file with the densities and the temepratures
   */
  void readInputFile();

  /**
   * Compute the initial solution for the coronal magnetic field using PFSS technique 
   */
  //void computePFSSMagneticField(const RealVector& stateCoordsSpherical,
                //RealVector& BPFSSCartesian,
                    //RealMatrix& sphCarTransMat);

  /**
   * Compute the Parker's solution for the solar wind 
   */
  //void computeParkerSolution(const CFreal r,
                 //RealVector& velParkerSpherical,
                 //CFreal& rhoParker,
                 //CFreal& pParker);
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
   * Construct the file path from which to read the data
   */
  boost::filesystem::path constructFilename(std::string fileName);

protected: // data

  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > _varSet;

  /// socket for the ions density
  Framework::DataSocketSource<std::vector<CFreal> > socket_IonsDens;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Alm, to be used in the boundary conditions
  //Framework::DataSocketSource<std::vector<CFreal> > socket_Almimg;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Blm, to be used in the boundary conditions
  //Framework::DataSocketSource<std::vector<CFreal> > socket_Blmreal;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Blm, to be used in the boundary conditions
  //Framework::DataSocketSource<std::vector<CFreal> > socket_Blmimg;

  /// name of the input file containing the data
  std::string _nameDataFile;

}; // class MFMHDInterpInitState2D2Fluid

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MFMHDInterpInitState2D2Fluid_hh

