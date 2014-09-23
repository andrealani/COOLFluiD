#ifndef COOLFluiD_Numerics_SpectralFV_SpectralFVMethodData_hh
#define COOLFluiD_Numerics_SpectralFV_SpectralFVMethodData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/FaceToCellGEBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/MeshData.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/NumericalJacobian.hh"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/VolumeIntegrator.hh"

#include "SpectralFV/CellToFaceGEBuilder.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

// forward declarations
class BCStateComputer;
class BaseBndFaceTermComputer;
class BaseFaceTermComputer;
class BaseVolTermComputer;
class FaceDiffusiveFlux;
class ReconstructStatesSpectralFV;
class RiemannFlux;
class SpectralFVElementData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents data object accessed by different SpectralFVMethodCom's
class SpectralFVMethodData : public Framework::SpaceMethodData {

public: // static functions

  /// Gets the Class name
  static std::string getClassName()
  {
    return "SpectralFVMethod";
  }

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

public: // functions

  /// Constructor
  SpectralFVMethodData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~SpectralFVMethodData();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the numerical jacobian calculator
   */
  Common::SafePtr<Framework::NumericalJacobian> getNumericalJacobian() const
  {
    cf_assert(m_numJacob.get() != CFNULL);
    return m_numJacob.get();
  }

  /**
   * Sets the ConvergenceMethod for this SpaceMethod to use
   * @pre the pointer to ConvergenceMethod is not constant to
   *      allow dynamic_casting
   */
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd)
  {
    m_convergenceMtd = convMtd;
  }

  /**
   * Get the ConvergenceMethod
   */
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> getConvergenceMethod() const
  {
    cf_assert(m_convergenceMtd.isNotNull());
    return m_convergenceMtd;
  }

  /// @return the GeometricEntity builder
  Common::SafePtr<
    Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > >
    getStdTrsGeoBuilder()
  {
    return &m_stdTrsGeoBuilder;
  }

  /// @return the GeometricEntity face builder
  Common::SafePtr<
    Framework::GeometricEntityPool< Framework::FaceToCellGEBuilder > >
    getFaceBuilder()
  {
    return &m_faceBuilder;
  }

  /// @return the GeometricEntity cell builder
  Common::SafePtr<
      Framework::GeometricEntityPool< SpectralFV::CellToFaceGEBuilder > >
      getCellBuilder()
  {
    return &m_cellBuilder;
  }

  /// @return the second GeometricEntity cell builder
  Common::SafePtr<
      Framework::GeometricEntityPool< SpectralFV::CellToFaceGEBuilder > >
      getSecondCellBuilder()
  {
    return &m_cellBuilder2nd;
  }

  /// @return the states reconstructor
  Common::SafePtr< ReconstructStatesSpectralFV > getStatesReconstructor();

  /// @return the boundary face term computer
  Common::SafePtr< BaseBndFaceTermComputer > getBndFaceTermComputer();

  /// @return an additional boundary face term computer
  Common::SafePtr< BaseBndFaceTermComputer > getAdditionalBndFaceTermComputer(const CFuint idx);

  /// @return the additional boundary face term computers
  Common::SafePtr< std::vector< Common::SafePtr< BaseBndFaceTermComputer > > > getAdditionalBndFaceTermComputers();

  /// @return the face term computer
  Common::SafePtr< BaseFaceTermComputer > getFaceTermComputer();

  /// @return an additional face term computer
  Common::SafePtr< BaseFaceTermComputer > getAdditionalFaceTermComputer(const CFuint idx);

  /// @return the additional face term computers
  Common::SafePtr< std::vector< Common::SafePtr< BaseFaceTermComputer > > > getAdditionalFaceTermComputers();

  /// @return the volume term computer
  Common::SafePtr< BaseVolTermComputer > getVolTermComputer();

  /// @return a second volume term computer
  Common::SafePtr< BaseVolTermComputer > getSecondVolTermComputer();

  /// @return the BCStateComputers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > getBCStateComputers();

  /// @return SafePtr to the Riemann flux strategy
  Common::SafePtr< RiemannFlux > getRiemannFlux();

  /// @return SafePtr to the FaceDiffusiveFlux flux strategy
  Common::SafePtr< FaceDiffusiveFlux > getFaceDiffusiveFlux();

  /// @return reference to m_bndFacesStartIdxs
  std::map< std::string , std::vector< std::vector< CFuint > > >& getBndFacesStartIdxs()
  {
    return m_bndFacesStartIdxs;
  }

  /// @return reference to m_innerFacesStartIdxs
  std::vector< CFuint >& getInnerFacesStartIdxs()
  {
    return m_innerFacesStartIdxs;
  }

  /// @return reference to m_svLocalData
  std::vector< SpectralFVElementData* >& getSVLocalData()
  {
    return m_svLocalData;
  }

  /// @return m_linearVarStr
  std::string getLinearVarStr()
  {
    return m_linearVarStr;
  }

  /// @return m_riemannFluxStr
  std::string getRiemannFluxStr()
  {
    return m_riemannFluxStr;
  }

  /// @return m_maxNbrStatesData
  CFuint getMaxNbrStatesData()
  {
    return m_maxNbrStatesData;
  }

  /// @return m_maxNbrStatesData
  CFuint getMaxNbrRFluxPnts()
  {
    return m_maxNbrRFluxPnts;
  }

  /// @return m_hasDiffTerm
  bool hasDiffTerm()
  {
    return m_hasDiffTerm;
  }

  /// set m_maxNbrStatesData
  void setMaxNbrStatesData(CFuint maxNbrStatesData)
  {
    m_maxNbrStatesData = maxNbrStatesData;
  }

  /// set m_maxNbrRFluxPnts
  void setMaxNbrRFluxPnts(CFuint maxNbrRFluxPnts)
  {
    m_maxNbrRFluxPnts = maxNbrRFluxPnts;
  }

  /// set m_resFactor
  void setResFactor(CFreal resFactor)
  {
    m_resFactor = resFactor;
  }

  /// @return m_resFactor
  CFreal getResFactor()
  {
    return m_resFactor;
  }

  /// @return m_createVolumesSocketBool
  bool createVolumesSocket()
  {
    return m_createVolumesSocketBool;
  }

  /// @return m_bcNameStr
  std::vector< std::string >& getBCNameStr()
  {
    return m_bcNameStr;
  }

  /// @return m_bcTRSNameStr
  Common::SafePtr< std::vector< std::vector< std::string > > > getBCTRSNameStr()
  {
    return &m_bcTRSNameStr;
  }

  /// Get the VolumeIntegrator
  Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator();

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

private:  // functions

  /**
   * Configures the Riemann flux
   */
  void configureRiemannFlux( Config::ConfigArgs& args );

  /**
   * Configures the face diffusive flux
   */
  void configureFaceDiffusiveFlux( Config::ConfigArgs& args );

  /**
   * Configures the boundary face term computer
   */
  void configureBndFaceTermComputer( Config::ConfigArgs& args );

  /**
   * Configures the face term computer
   */
  void configureFaceTermComputer( Config::ConfigArgs& args );

  /**
   * Configures the volume term computer
   */
  void configureVolTermComputer( Config::ConfigArgs& args );

  /**
   * Configures the BC state computers
   */
  void configureBCStateComputers( Config::ConfigArgs& args );

  /**
   * Creates the local data for SV
   */
  void createSVLocalData();

private:  // data

  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> m_numJacob;

  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> m_convergenceMtd;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder >  m_stdTrsGeoBuilder;

  /// Builder for faces (containing the neighbouring cells)
  Framework::GeometricEntityPool< Framework::FaceToCellGEBuilder >  m_faceBuilder;

  /// Builder for cells (containing the neighbouring faces)
  Framework::GeometricEntityPool< SpectralFV::CellToFaceGEBuilder >  m_cellBuilder;

  /// Second builder for cells (containing the neighbouring faces)
  Framework::GeometricEntityPool< SpectralFV::CellToFaceGEBuilder >  m_cellBuilder2nd;

  /// pointer to states reconstructor strategy
  Common::SelfRegistPtr< ReconstructStatesSpectralFV > m_statesReconstructor;

  /// String for the volume terms computer
  std::string m_volTermComputerStr;

  /// pointer to volume term computer
  Common::SelfRegistPtr< BaseVolTermComputer > m_volTermComputer;

  /// pointer to second volume term computer
  Common::SelfRegistPtr< BaseVolTermComputer > m_volTermComputer2nd;

  /// String for the face terms computer
  std::string m_faceTermComputerStr;

  /// pointer to face term computer
  Common::SelfRegistPtr< BaseFaceTermComputer > m_faceTermComputer;

  /// pointers to additional face term computers
  std::vector< Common::SelfRegistPtr< BaseFaceTermComputer > > m_addFaceTermComputers;

  /// SafePointers to additional face term computers
  std::vector< Common::SafePtr< BaseFaceTermComputer > > m_addFaceTermComputersSP;

  /// String for the boundary face terms computer
  std::string m_bndFaceTermComputerStr;

  /// pointer to boundary face term computer
  Common::SelfRegistPtr< BaseBndFaceTermComputer > m_bndFaceTermComputer;

  /// pointers to additional boundary face term computer
  std::vector< Common::SelfRegistPtr< BaseBndFaceTermComputer > > m_addBndFaceTermComputers;

  /// SafePointers to additional boundary face term computer
  std::vector< Common::SafePtr< BaseBndFaceTermComputer > > m_addBndFaceTermComputersSP;

  /// String for the linear variable name (for instance, the Roe average variables)
  std::string m_linearVarStr;

  /// String for the Riemann flux
  std::string m_riemannFluxStr;

  /// String for the face diffusive flux
  std::string m_faceDiffFluxStr;

  /// pointer to Riemann flux strategy
  Common::SelfRegistPtr< RiemannFlux > m_riemannFlux;

  /// pointer to FaceDiffusiveFlux strategy
  Common::SelfRegistPtr< FaceDiffusiveFlux > m_faceDiffFlux;

  /// The boundary condition state computer strategies
  std::vector< Common::SelfRegistPtr< BCStateComputer > > m_bcs;

  /// The boundary condition state computer strategies, as SafePtrs
  std::vector< Common::SafePtr< BCStateComputer > > m_bcsSP;

  /// The boundary condition strategy types
  std::vector< std::string > m_bcTypeStr;

  /// The boundary condition strategy names for configuration
  std::vector< std::string > m_bcNameStr;

  /// The boundary condition TRS names
  std::vector< std::vector< std::string > > m_bcTRSNameStr;

  /// start index of inner faces with a certain orientation
  std::vector< CFuint > m_innerFacesStartIdxs;

  /// map between the boundary TRS and the start index of faces with a certain orientation
  std::map< std::string , std::vector< std::vector< CFuint > > > m_bndFacesStartIdxs;

  /// vector containing the  SpectralFVElementData for different element types
  std::vector< SpectralFVElementData* > m_svLocalData;

  /// variable for maximum number of statesData that has to be computed
  CFuint m_maxNbrStatesData;

  /// variable for maximum number of points in which the Riemann solver is evaluated
  CFuint m_maxNbrRFluxPnts;

  /// boolean storing wether there is a diffusive term
  bool m_hasDiffTerm;

  /// factor to multiply the residual with, coming from the time discretization
  CFreal m_resFactor;

  /// boolean telling wheter the socket containing the volume for each state (!= cell) has to be created
  bool m_createVolumesSocketBool;

};  // end of class SpectralFVMethodData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a MethodCommand for SpectralFV
typedef Framework::MethodCommand< SpectralFVMethodData > SpectralFVMethodCom;

/// Definition of a command provider for SpectralFV
typedef SpectralFVMethodCom::PROVIDER SpectralFVMethodComProvider;

/// Definition of a MethodStrategy for SpectralFV
typedef Framework::MethodStrategy< SpectralFVMethodData > SpectralFVMethodStrategy;

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SpectralFV_SpectralFVMethodData_hh

