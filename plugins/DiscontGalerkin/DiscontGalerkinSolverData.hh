
#ifndef COOLFluiD_DiscontGalerkin_DiscontGalerkinSolverData_hh
#define COOLFluiD_DiscontGalerkin_DiscontGalerkinSolverData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntityPool.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/ContourIntegrator.hh"

#include "Framework/DofDataHandleIterator.hh"
#include "Framework/ProxyDofIterator.hh"

#include "Framework/FaceToCellGEBuilder.hh"

#include "Common/StringOps.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//     class DiscontGalerkinStrategy;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents data object accessed by different DiscontGalerkinSolverCom's
 * @author Tiago Quintino
 * @author Pedro Maciel
 */
class DiscontGalerkinSolverData : public Framework::SpaceMethodData {

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  DiscontGalerkinSolverData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~DiscontGalerkinSolverData();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

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
    return &m_FaceBuilder;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "DiscontGalerkinSolver";
  }

  /// Get the VolumeIntegrator
  Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator();

  /// Get the ContourIntegrator
  Common::SafePtr< Framework::ContourIntegrator > getContourIntegrator();

  /// Sets up the DiscontGalerkinData
  void setup();

  /// Gets the empty strategy
//   DiscontGalerkinStrategy * getDiscontGalerkinStrategy() const
//   {
//     cf_assert(m_emptystrategy.isNotNull());
//     return m_emptystrategy.getPtr();
//   }

  /// Get the max eigenvalue
  CFreal getMaxEigenval()
  {
    return m_maxEigenValue;
  }

  /// Get the max eigenvalue
  void setMaxEigenval(CFreal eigenval)
  {
    m_maxEigenValue=eigenval;
  }

  /// Get sigma constant
  CFreal getSigma()
  {
    return(Common::StringOps::from_str < CFreal > ( m_sigmaStr ));
  }

  /// Get theta constant
  CFreal getTheta()
  {
    if (m_typeStr == "IIPG")
      return 0;
    else if (m_typeStr == "NIPG")
      return -1;
    else if (m_typeStr == "SIPG")
      return 1;
    else 
    {
      CFout << "ERROR type of DGFEM not set. Using IIPG." << CFendl;
      return 0;
    }
  }

  /// Get reynolds number
  CFreal getReynolds()
  {
    cf_assert(m_reynoldStr != "");
    return Common::StringOps::from_str < CFreal > ( m_reynoldStr );
  }

  /// Get max CFL constant
  CFreal getMaxCFL()
  {
    return Common::StringOps::from_str < CFreal > ( m_maxCFLStr );
  }

  /// Get alpha constant
  CFreal getAlpha()
  {
    return Common::StringOps::from_str < CFreal > ( m_alphaStr );
  }

  /// Get residual
  CFreal getResidual()
  {
    return(m_Residual);
  }
  /// Set residual
void setResidual(CFreal res)
  {
    m_Residual = res;
    return;
  }

private:  // helper functions

  /// Configures the ContourIntegrator and the IntegrableEntity
  void configureIntegrator();

private:  // data

  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> m_convergenceMtd;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > m_stdTrsGeoBuilder;

  /// The volume integrator
  Framework::VolumeIntegrator m_volumeIntegrator;

  /// The volume integrator
  Framework::ContourIntegrator m_contourIntegrator;

  /// String for configuring the numerical integrator QuadratureType
  std::string m_volintquadStr;

  /// String for configuring the numerical integrator QuadratureType
  std::string m_conintquadStr;

  /// String for configuring the numerical integrator Order
  std::string m_volintorderStr;

  /// String for `configuring the numerical integrator Order
  std::string m_conintorderStr;

  /// String for configuring the type of the DGFEM - IIPG(0), SIPG(1), NIPG(-1)
  std::string m_typeStr;

  /// String for configuring the constant reynolds in viscous flow
  std::string m_reynoldStr;

  /// String for setting max of CFL for computing of time step
  std::string m_maxCFLStr;

  /// String for configuring the constant alpha for computing of time step
  std::string m_alphaStr;

  /// String for configuring the constant sigma in viscous flow
  std::string m_sigmaStr;

  /// Builder for spectral faces
  Framework::GeometricEntityPool< Framework::FaceToCellGEBuilder >  m_FaceBuilder;

  ///store of max eigenvalue
  CFreal m_maxEigenValue;

  ///store of residual
  CFreal m_Residual;

  /// DiscontGalerkin strategy
//   Common::SelfRegistPtr< DiscontGalerkinStrategy > m_emptystrategy;

  /// String to configure empty strategy
//   std::string m_emptystrategyStr;

};  // end of class DiscontGalerkinSolverData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a MethodCommand for DiscontGalerkin
typedef Framework::MethodCommand< DiscontGalerkinSolverData > DiscontGalerkinSolverCom;

/// Definition of a MethodStrategy for DiscontGalerkin
typedef Framework::MethodStrategy< DiscontGalerkinSolverData > DiscontGalerkinSolverStrategy;

/// Definition of a command provider for DiscontGalerkin
typedef Framework::MethodCommand<DiscontGalerkinSolverData>::PROVIDER DiscontGalerkinSolverComProvider;

//////////////////////////////////////////////////////////////////////////////

  } // namespace DiscontGalerkin
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_DiscontGalerkinSolverData_hh

