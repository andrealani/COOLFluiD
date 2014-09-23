#ifndef COOLFluiD_Numerics_FiniteElement_FiniteElementMethodData_hh
#define COOLFluiD_Numerics_FiniteElement_FiniteElementMethodData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"

#include "Framework/MethodCommand.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ComputeCFL.hh"
#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/InertiaVarSet.hh"
#include "Framework/SourceVarSet.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/SpaceMethodData.hh"
#include "FEM_VolumeIntegrator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

    class DiffusiveEntity;
    class InertiaEntity;
    class ConvectiveEntity;
    class LinearSourceEntity;
    class IndepSourceEntity;
    class ComputeJacobStrategy;
    class ComputeResidualStrategy;
    class ComputeConvectiveTerm;
    class ComputeDiffusiveTerm;
    class ComputeLinearSourceTerm;
    class ComputeIndepSourceTerm;
    class ComputeInertiaTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Data Object that is accessed by the different
 * FiniteElementMethodCom 's that compose the FiniteElementMethod.
 *
 * @todo there is missing documentation in this class.
 *
 * @author Tiago Quintino
 *
 */
class FiniteElementMethodData : public Framework::SpaceMethodData {

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  FiniteElementMethodData(Common::SafePtr<Framework::Method> owner);

  /**
   * Destructor.
   */
  ~FiniteElementMethodData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /**
   * Get the inertia variables set
   */
  Common::SafePtr<Framework::InertiaVarSet> getInertiaVar() const
  {
    cf_assert(_inertiaVar.isNotNull());
    return _inertiaVar.getPtr();
  }

  /**
   * Get the linear Source variables set
   */
  Common::SafePtr<Framework::SourceVarSet> getSourceVar() const
  {
    cf_assert(_sourceVar.isNotNull());
    return _sourceVar.getPtr();
  }

  /**
   * Sets the ConvergenceMethod for this SpaceMethod to use
   * @pre the pointer to ConvergenceMethod is not constant to
   *      allow dynamic_casting
   */
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd)
  {
    _convergenceMtd = convMtd;
  }

  /**
   * Get the ConvergenceMethod
   */
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> getConvergenceMethod() const
  {
    cf_assert(_convergenceMtd.isNotNull());
    return _convergenceMtd;
  }

  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
  getStdTrsGeoBuilder()
  {
    return &_stdTrsGeoBuilder;
  }

  /**
   * Set the residual factor
   */
  void setResFactor(const CFreal& factor)
  {
    _resFactor = factor;
  }

  /**
   * Get the residual factor
   */
  CFreal getResFactor() const
  {
    return _resFactor;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "FiniteElementMethod";
  }

  /**
   * Gets the convective term computer
   */
  Common::SafePtr<ComputeConvectiveTerm> getConvectiveTermComputer()
  {
    return _convectiveTerm.getPtr();
  }

  /**
   * Gets the diffusive term computer
   */
  Common::SafePtr<ComputeDiffusiveTerm> getDiffusiveTermComputer()
  {
    return _diffusiveTerm.getPtr();
  }

  /**
   * Gets the diffusive entity
   */
  Common::SafePtr<DiffusiveEntity> getDiffusiveEntity()
  {
    return _diffusiveEntity.getPtr();
  }

  /**
   * Gets the inertia entity
   */
  Common::SafePtr<InertiaEntity> getInertiaEntity()
  {
    return _inertiaEntity.getPtr();
  }

  /**
   * Gets the diffusive entity
   */
  Common::SafePtr<ConvectiveEntity> getConvectiveEntity()
  {
    return _convectiveEntity.getPtr();
  }

  /**
   * Gets the independent source term computer
   */
  Common::SafePtr<ComputeIndepSourceTerm> getIndepSourceTermComputer()
  {
    return _indepSourceTerm.getPtr();
  }

  /**
   * Gets the Independent source entity
   */
  Common::SafePtr<IndepSourceEntity> getIndepSourceEntity()
  {
    return _indepSourceEntity.getPtr();
  }

  /**
   * Gets the linear source term computer
   */
  Common::SafePtr<ComputeLinearSourceTerm> getLinearSourceTermComputer()
  {
    return _linearSourceTerm.getPtr();
  }

  /**
   * Gets the Independent source entity
   */
  Common::SafePtr<LinearSourceEntity> getLinearSourceEntity()
  {
    return _linearSourceEntity.getPtr();
  }

  /**
   * Gets the inertia term computer
   */
  Common::SafePtr<ComputeInertiaTerm> getInertiaTermComputer()
  {
    return _inertiaTerm.getPtr();
  }

  /**
   * Gets  the strategy to compute the system jacobian matrix
   */
  Common::SafePtr<ComputeJacobStrategy> getJacobianStrategy() const
  {
    cf_assert(_jacobianStrategy.isNotNull());
    return _jacobianStrategy.getPtr();
  }

  /**
   * Gets  the strategy to compute the element residual
   */
  Common::SafePtr<ComputeResidualStrategy> getResidualStrategy() const
  {
    cf_assert(_residualStrategy.isNotNull());
    return _residualStrategy.getPtr();
  }

  /**
   * Get the VolumeIntegrator
   */
  Common::SafePtr<Framework::VolumeIntegrator> getVolumeIntegrator();

  /**
   * Get the FEM_VolumeIntegrator
   */
  Common::SafePtr<FEM_VolumeIntegrator> getFEMVolumeIntegrator();

  /**
   * Check if a Dirichlet BC has already been applied (you should not apply a NeumannBC after a DirichletBC)
   */
  bool isDirichletBCApplied() const
  {
    return _isDirichletBCApplied;
  }

  /**
   * Set flag to know if a Dirichlet BC has already been applied
   * (you should not apply a NeumannBC after a DirichletBC)
   */
  void setDirichletBCApplied(bool applied)
  {
    _isDirichletBCApplied = applied;
  }

  /**
   * Get the Local Element Data
   */
  LocalElementData& getLocalElementData()
  {
    return m_local_elem_data;
  }

private: // functions

  /**
   * Configures the ContourIntegrator and the IntegrableEntity
   */
  void configureIntegrator ( Config::ConfigArgs& args );

  /**
   * Configures the term computers
   */
  void configureTerms ( Config::ConfigArgs& args );

  /**
   * Configures the VarSet's
   */
  void configureVarSets ( Config::ConfigArgs& args );

  /**
   * Configures the variable set transfomers
   */
  void configureVarSetTransformers ( Config::ConfigArgs& args );

private: // data

  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> _convergenceMtd;

  /// Inertia variables set
  Common::SelfRegistPtr<Framework::InertiaVarSet> _inertiaVar;

  /// Linear Source variables set
  Common::SelfRegistPtr<Framework::SourceVarSet> _sourceVar;

  /// Convective term computer
  Common::SelfRegistPtr<ComputeConvectiveTerm> _convectiveTerm;

  /// Convective Entity
  Common::SelfRegistPtr<ConvectiveEntity> _convectiveEntity;

  /// Diffusive term computer
  Common::SelfRegistPtr<ComputeDiffusiveTerm> _diffusiveTerm;

  /// Diffusive Entity
  Common::SelfRegistPtr<DiffusiveEntity> _diffusiveEntity;

  /// Inertia term computer
  Common::SelfRegistPtr<ComputeInertiaTerm> _inertiaTerm;

  /// Inertia Entity
  Common::SelfRegistPtr<InertiaEntity> _inertiaEntity;

  /// Linear Source term computer
  Common::SelfRegistPtr<ComputeLinearSourceTerm> _linearSourceTerm;

  /// Linear Source Entity
  Common::SelfRegistPtr<LinearSourceEntity> _linearSourceEntity;

  /// Independent Source term computer
  Common::SelfRegistPtr<ComputeIndepSourceTerm> _indepSourceTerm;

  /// Independent Source Entity
  Common::SelfRegistPtr<IndepSourceEntity> _indepSourceEntity;

  /// strategy to compute the jacobian system matrix
  Common::SelfRegistPtr<ComputeJacobStrategy> _jacobianStrategy;

  /// strategy to compute the residual
  Common::SelfRegistPtr<ComputeResidualStrategy> _residualStrategy;

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;

  /// residual factor
  CFreal _resFactor;

  /// string for the configuration of the numerical integrator QuadratureType
  std::string _integratorQuadratureStr;

  /// string for the configuration of the numerical integrator Order
  std::string _integratorOrderStr;

  /// inertia variable set name
  std::string _inertiaVarStr;

  /// source variable set name
  std::string _sourceVarStr;

  /// independent source entity name
  std::string  _indepSourceEntityStr;

  /// string to configure _jacobianStrategy
  std::string _jacobianStrategyStr;

  /// string to configure _residualStrategy
  std::string _residualStrategyStr;

  /// The FEM volume Integrator
  FEM_VolumeIntegrator _femVolumeIntegrator;

  ///Flag to know if a dirichlet BC has already been applied
  bool _isDirichletBCApplied;

  /// Data relative to the current element being processed
  LocalElementData m_local_elem_data;

}; // end of class FiniteElementMethodData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a MethodCommand for FiniteElement
typedef Framework::MethodCommand<FiniteElementMethodData> FiniteElementMethodCom;

/// Definition of a MethodCommand provider for FiniteElement
typedef Framework::MethodCommand<FiniteElementMethodData>::PROVIDER FiniteElementMethodComProvider;

/// Definition of a Strategy for FiniteElement
typedef Framework::MethodStrategy<FiniteElementMethodData> FiniteElementMethodStrategy;

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_FiniteElementMethodData_hh
