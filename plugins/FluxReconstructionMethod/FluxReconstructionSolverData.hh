// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_FluxReconstructionSolverData_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_FluxReconstructionSolverData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntityPool.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/VolumeIntegrator.hh"

#include "Framework/DofDataHandleIterator.hh"
#include "Framework/ProxyDofIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

    class FluxReconstructionStrategy;
    class BaseInterfaceFlux;

//////////////////////////////////////////////////////////////////////////////

/// This class represents data object accessed by different FluxReconstructionSolverCom's
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class FluxReconstructionSolverData : public Framework::SpaceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  FluxReconstructionSolverData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~FluxReconstructionSolverData();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre LinearSystemSolver pointer is not constant to allow dynamic_casting
  void setLinearSystemSolver(
    Framework::MultiMethodHandle< Framework::LinearSystemSolver > lss )
  {
    m_lss = lss;
  }

  /// Get the linear system solver
  Framework::MultiMethodHandle< Framework::LinearSystemSolver >
    getLinearSystemSolver() const
  {
    cf_assert(m_lss.isNotNull());
    return m_lss;
  }

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to
  ///      allow dynamic_casting
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd)
  {
    m_convergenceMtd = convMtd;
  }

  /// Get the ConvergenceMethod
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

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FluxReconstructionSolver";
  }

  /// Get the VolumeIntegrator
  Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator();

  /// Sets up the FluxReconstructionData
  void setup();

  /// Gets the FluxReconstruction strategy
  FluxReconstructionStrategy * getFluxReconstructionStrategy() const
  {
    cf_assert(m_fluxreconstructionstrategy.isNotNull());
    return m_fluxreconstructionstrategy.getPtr();
  }
  
  /// Gets the interface flux computation strategy
  BaseInterfaceFlux * getInterfaceFluxStrategy() const
  {
    cf_assert(m_interfaceflux.isNotNull());
    return m_interfaceflux.getPtr();
  }

private:  // helper functions

  /// Configures the ContourIntegrator and the IntegrableEntity
  void configureIntegrator();

private:  // data

  /// Linear system solver
  Framework::MultiMethodHandle< Framework::LinearSystemSolver > m_lss;

  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> m_convergenceMtd;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > m_stdTrsGeoBuilder;

  /// The volume integrator
  Framework::VolumeIntegrator m_volumeIntegrator;

  /// String for configuring the numerical integrator QuadratureType
  std::string m_intquadStr;

  /// String for configuring the numerical integrator Order
  std::string m_intorderStr;

  /// FluxReconstruction strategy
  Common::SelfRegistPtr< FluxReconstructionStrategy > m_fluxreconstructionstrategy;

  /// String to configure FluxReconstruction strategy
  std::string m_fluxreconstructionstrategyStr;
  
  /// Interface flux computation strategy
  Common::SelfRegistPtr< BaseInterfaceFlux > m_interfaceflux;

  /// String to configure interface flux computation strategy
  std::string m_interfacefluxStr;

};  // end of class FluxReconstructionSolverData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a MethodCommand for FluxReconstructionMethod
typedef Framework::MethodCommand< FluxReconstructionSolverData > FluxReconstructionSolverCom;

/// Definition of a command provider for FluxReconstructionMethod
typedef FluxReconstructionSolverCom::PROVIDER FluxReconstructionSolverComProvider;

/// Definition of a MethodStrategy for FluxReconstructionMethod
typedef Framework::MethodStrategy< FluxReconstructionSolverData > FluxReconstructionSolverStrategy;

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_FluxReconstructionSolverData_hh

