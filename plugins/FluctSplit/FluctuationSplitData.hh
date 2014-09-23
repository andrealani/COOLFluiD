#ifndef COOLFluiD_Numerics_FluctSplit_FluctuationSplitData_hh
#define COOLFluiD_Numerics_FluctSplit_FluctuationSplitData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Common/SafePtr.hh"
#include "Config/ConfigObject.hh"

#include "Framework/MethodCommand.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/ComputeCFL.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/VolumeIntegrator.hh"

#include "FluctSplit/SourceTermSplitter.hh"
#include "FluctSplit/DistributionData.hh"
#include "FluctSplit/ComputeJacobianFix.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

      class FluctuationSplitStrategy;
      class ComputeJacobStrategy;
      class ComputeDiffusiveTerm;
      class ArtificialDiffusionStrategy;
      class ComputeSourceTermFSM;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// FluctuationSplitCom 's that compose the FluctuationSplit.
/// @todo there is missing documentation in this class.
/// @author Tiago Quintino
/// @author Andrea Lani
class FluctSplit_API FluctuationSplitData : public Framework::SpaceMethodData {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer VectorTransformer;

  /// Default constructor without arguments.
  FluctuationSplitData(Common::SafePtr<Framework::Method> owner);

  /// Destructor.
  ~FluctuationSplitData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Set the residual factor
  void setResFactor(const CFreal& factor) { m_resFactor = factor; }

  /// Get the residual factor
  CFreal getResFactor() const {  return m_resFactor; }

  /// @return true if using both System and Scalar splitter
  bool isMultipleSplitter() const
  {
    cf_assert(m_SysSplitter.isNotNull());
    cf_assert(m_SclSplitter.isNotNull());
    return m_multipleSplitter;
  }

  /// @return true if using only one splitter either system or scalar.
  bool isSingleSplitter() const
  {
    return !isMultipleSplitter();
  }

  /// @return true if using only System splitter
  bool isOnlySysSplitter() const
  {
    cf_assert(m_SysSplitter.isNotNull());
    return m_isOnlySysSplitter;
  }

  /// @return true if the SubSystem is axisymmetric
  bool isAxisymmetric() const
  {
    return m_isAxisymm;
  }

  /// @return true if the source term has to be added to the convective flux
  bool includeSourceInFlux() const
  {
    return m_includeSourceInFlux;
  }

  /// Get System Splitter
  /// @return the pointer to the system splitter
  Common::SafePtr<Splitter> getSysSplitter() const {  return m_SysSplitter;  }

  /// Get Scalar Splitter
  /// @return the pointer to the system splitter
  Common::SafePtr<Splitter> getScalarSplitter() const  {  return m_SclSplitter;  }

  /// Get the system source term Splitter
  /// @return the pointer to the system source term splitter
  Common::SafePtr<std::vector<Common::SelfRegistPtr<SourceTermSplitter> > > getSourceTermSplitter()
  {
    cf_assert(m_sourceTermSplitter.size() > 0);
    return &m_sourceTermSplitter;
  }

  /// Get the system source term Splitter
  /// @return the pointer to the system source term splitter
  Common::SafePtr<SourceTermSplitter> getSourceTermSplitter(CFuint is)
  {
    cf_assert(is < m_sourceTermSplitter.size());
    return m_sourceTermSplitter[is].getPtr();
  }

  /// Get Splitter
  /// @return the pointer to the system splitter
  Common::SafePtr<Splitter> getSplitter() const
  {
    cf_assert(!isMultipleSplitter());

    if (m_SclSplitter->isNull()) {
      // cf_assert(!m_SysSplitter->isNull());
      return m_SysSplitter;
    }
    return m_SclSplitter;
  }

  /// Switch the splitters to computation of implicit jacobian
  void switchToImplicit()
  {
    m_SysSplitter = m_jacobSysSplitter.getPtr();
    m_SclSplitter = m_jacobSclSplitter.getPtr();
  }

  /// Switch the splitters to explicit computation or RHS computation
  void switchToExplicit()
  {
    m_SysSplitter = m_rhsSysSplitter.getPtr();
    m_SclSplitter = m_rhsSclSplitter.getPtr();
  }

  /// Get the diffusive term computer
  Common::SafePtr<ComputeDiffusiveTerm> getDiffusiveTermComputer() const
  {
    cf_assert(m_diffTermComputer.isNotNull());
    return m_diffTermComputer.getPtr();
  }

  /// Get the distribution variables set
  Common::SafePtr<Framework::ConvectiveVarSet> getDistribVar() const
  {
    cf_assert(m_distribVar.isNotNull());
    return m_distribVar.getPtr();
  }

  /// Get the linearizion variables set
  Common::SafePtr<Framework::ConvectiveVarSet> getLinearVar() const
  {
    cf_assert(m_linearVar.isNotNull());
    return m_linearVar.getPtr();
  }

  /// Get Linearizer
  /// @return the pointer to the jacobian linearizer
  Common::SafePtr<Framework::JacobianLinearizer> getLinearizer() const
  {
    cf_assert(m_linearizer.isNotNull());
    return m_linearizer.getPtr();
  }

  /// Get the source term computer
  /// @return the pointer to the source term computer
  Common::SafePtr< std::vector<Common::SelfRegistPtr<ComputeSourceTermFSM> > >  getSourceTermComputer();

  /// Get the system source term Splitter
  /// @return the pointer to the system source term splitter
  Common::SafePtr<ComputeSourceTermFSM> getSourceTermComputer(CFuint is);

  /// Get the jacobian fix computer
  Common::SafePtr<FluctSplit::ComputeJacobianFix> getJacobianFixComputer() const
  {
    cf_assert(m_jacobFixComputer.isNotNull());
    return m_jacobFixComputer.getPtr();
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
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
  getStdTrsGeoBuilder()
  {
    return &m_stdTrsGeoBuilder;
  }

  /// Get the matrix transformation from Solution to Distribution variables
  /// @return the pointer to the transformation
  Common::SafePtr<MatrixTransformer> getSolutionToDistribMatTrans() const
  {
    cf_assert(m_solutionToDistMatTrans.isNotNull());
    return m_solutionToDistMatTrans.getPtr();
  }

  /// Get the matrix transformation from Distribution to Solution variables
  /// @return the pointer to the transformation
  Common::SafePtr<MatrixTransformer> getDistribToSolutionMatTrans() const
  {
    cf_assert(m_distToSolutionMatTrans.isNotNull());
    return m_distToSolutionMatTrans.getPtr();
  }

  /// Get the matrix transformation from Linear to Distribution variables
  /// @return the pointer to the transformation
  Common::SafePtr<MatrixTransformer> getLinearToDistribMatTrans() const
  {
    cf_assert(m_linearToDistMatTrans.isNotNull());
    return m_linearToDistMatTrans.getPtr();
  }

  /// Get the matrix transformation from Solution to Linear in Update
  /// @return the pointer to the transformation
  Common::SafePtr<MatrixTransformer> getSolutionToLinearInUpdateMatTrans() const
  {
    cf_assert(m_solutionToLinearInUpdateMatTrans.isNotNull());
    return m_solutionToLinearInUpdateMatTrans.getPtr();
  }

  /// Get the matrix transformation from Solution to Linear variables
  /// @return the pointer to the transformation
  Common::SafePtr<MatrixTransformer> getSolutionToLinearMatTrans() const
  {
    cf_assert(m_solutionToLinearMatTrans.isNotNull());
    return m_solutionToLinearMatTrans.getPtr();
  }

  /// Get the vector transformation from Update to Linear variables
  /// @return the pointer to the transformation
  Common::SafePtr<VectorTransformer> getUpdateToLinearVecTrans() const
  {
    cf_assert(m_updateToLinearVecTrans.isNotNull());
    return m_updateToLinearVecTrans.getPtr();
  }

  /// Get the vector transformation from Update to Solution variables
  /// @return the pointer to the transformation
  Common::SafePtr<VectorTransformer> getUpdateToSolutionVecTrans() const
  {
    cf_assert(m_updateToSolutionVecTrans.isNotNull());
    return m_updateToSolutionVecTrans.getPtr();
  }

  /// Get the matrix transformer from solution to update variables
  /// starting from update variables
  Common::SafePtr<Framework::VarSetMatrixTransformer>
  getSolToUpdateInUpdateMatTrans() const
  {
    cf_assert(m_solToUpdateInUpdateMatTrans.isNotNull());
    return m_solToUpdateInUpdateMatTrans.getPtr();
  }

  /// Get the matrix transformer from update to solution variables
  /// starting from update variables
  Common::SafePtr<Framework::VarSetMatrixTransformer>
  getUpdateToSolutionInUpdateMatTrans() const
  {
    cf_assert(m_updateToSolutionInUpdateMatTrans.isNotNull());
    return m_updateToSolutionInUpdateMatTrans.getPtr();
  }

  /// Tell if a varaible has to be applied for the residual
  bool isResidualTransformationNeeded() const
  {
    return (_updateVarStr != _solutionVarStr);
  }

  /// Fluctuation splitting strategy name
  std::string getFluctSplitStrategyName() const
  {
    if (isMultipleSplitter()) {
      return m_fsStrategyName + "Multi";
    }
    return m_fsStrategyName;
  }

  /// Fluctuation splitting strategy name
  std::string getArtDiffStrategyName() const
  {

    return m_adStrategyName;
  }

  /// Get the Fluctuation splitting strategy
  Common::SafePtr<FluctuationSplitStrategy> getFluctSplitStrategy() const;

  /// Get the Fluctuation splitting strategy
  Common::SafePtr<ArtificialDiffusionStrategy> getArtificialDiffusionStrategy() const;

  /// Get the ComputeJacobStrategy
  Common::SafePtr<ComputeJacobStrategy> getJacobStrategy() const;

  /// Get the ContourIntegrator
  Common::SafePtr<Framework::ContourIntegrator> getContourIntegrator();

  /// Get the VolumeIntegrator
  Common::SafePtr<Framework::VolumeIntegrator> getVolumeIntegrator();

  /// Get the data for the fluctuation distribution
  DistributionData& getDistributionData() {  return m_distData;  }

  /// Get block separator index
  /// @return integer with block separator index
  CFuint getBlockSeparator() const;

  /// Jacobian computation strategy name
  std::string getComputeJacobianStrategyName() const
  {
    return m_jacobStrategyName;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FluctuationSplit";
  }
  
  /// Apply first the scalar and the n the system scheme (by default is the opposite)
  bool isScalarFirst() const
  {
    return m_scalarFirst;
  }
  
  /// Check is artificial diffusion term should be added
  bool hasArtificialDiff() const
  {
    return m_artdiff;
  }

  /// Set a flag to know if we are in the initialization phase
  void setInitializationPhase(bool flag)
  {
    m_isInitializationPhase = flag;

  }

  /// Check if we are currently in the initialization phase
  bool isInitializationPhase() const
  {
    return m_isInitializationPhase;
  }

private:

  /// Configures the VarSet's
  void configureVarSets ( Config::ConfigArgs& args );

  /// Configures the Transformers
  void configureTransformers ( Config::ConfigArgs& args );

  /// Configures the Linearizer
  void configureLinearizer ( Config::ConfigArgs& args );

  /// Configures the FluctuationSplitStrategy
  void configureFluctSplitStrategy ( Config::ConfigArgs& args );

  /// Configures the FluctuationSplitStrategy
  void configureArtDiffStrategy ( Config::ConfigArgs& args );

  /// Configures the FluctuationSplitStrategy
  void configureJacobStrategy ( Config::ConfigArgs& args );

  /// Configures the source term computer
  void configureSourceTermComputer ( Config::ConfigArgs& args );

  /// Configures the diffusive term computer
  void configureDiffusiveTermComputer ( Config::ConfigArgs& args );

  /// Configures the CFL computer
  void configureCFLComputer ( Config::ConfigArgs& args );

  /// Configures the Splitters
  /// @pre assumes you already called configureTransformers()
  void configureSplitters ( Config::ConfigArgs& args );

  /// Configures the Source Term Splitters
  /// @pre assumes you already called configureTransformers()
  void configureSourceTermSplitters ( Config::ConfigArgs& args );

  /// Configures both Integrators
  void configureIntegrators ( Config::ConfigArgs& args );

private:

  /// Multiple Splitter
  bool m_multipleSplitter;

  /// The splitter is only a system scheme
  bool m_isOnlySysSplitter;

  /// System Splitter for RHS and explicit computations
  Common::SafePtr<Splitter> m_SysSplitter;

  /// Scalar Splitter for RHS and explicit computations
  Common::SafePtr<Splitter> m_SclSplitter;

  /// System Splitter for RHS and explicit computations
  Common::SelfRegistPtr<Splitter> m_rhsSysSplitter;

  /// Scalar Splitter for RHS and explicit computations
  Common::SelfRegistPtr<Splitter> m_rhsSclSplitter;

  /// System Splitter for the Jacobian of the implicit computations
  Common::SelfRegistPtr<Splitter> m_jacobSysSplitter;

  /// Scalar Splitter for the Jacobian of the implicit computations
  Common::SelfRegistPtr<Splitter> m_jacobSclSplitter;

  /// system Source Term splitter
  std::vector<Common::SelfRegistPtr<SourceTermSplitter> > m_sourceTermSplitter;

  /// variable set in which distributing the residual
  Common::SelfRegistPtr<Framework::ConvectiveVarSet>  m_distribVar;

  /// Linearizer
  Common::SelfRegistPtr<Framework::JacobianLinearizer> m_linearizer;

  /// Source term computer
  std::vector<Common::SelfRegistPtr<ComputeSourceTermFSM> > m_stComputer;

  /// Diffusive term computer
  Common::SelfRegistPtr<ComputeDiffusiveTerm> m_diffTermComputer;

  /// Jacobian fix computer
  Common::SelfRegistPtr<ComputeJacobianFix> m_jacobFixComputer;

  /// Transformer from Distribution to Solution Variables
  Common::SelfRegistPtr<MatrixTransformer> m_solutionToDistMatTrans;

  /// Transformer from Solution to Distribution Variables
  Common::SelfRegistPtr<MatrixTransformer> m_distToSolutionMatTrans;

  /// Transformer from Solution to Distribution Variables
  Common::SelfRegistPtr<MatrixTransformer> m_linearToDistMatTrans;

  /// Transformer from Linear to Solution in Solution Variables
  Common::SelfRegistPtr<MatrixTransformer> m_solutionToLinearInUpdateMatTrans;

  /// Transformer from Solution to Linear Variables
  Common::SelfRegistPtr<MatrixTransformer> m_solutionToLinearMatTrans;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<VectorTransformer> m_updateToLinearVecTrans;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<VectorTransformer> m_updateToSolutionVecTrans;

  /// Matrix transformer from solution to update variables
  /// starting from update variables
  Common::SelfRegistPtr<MatrixTransformer> m_solToUpdateInUpdateMatTrans;

  /// Matrix transformer from update to solution variables
  /// starting from update variables
  Common::SelfRegistPtr<MatrixTransformer> m_updateToSolutionInUpdateMatTrans;

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> m_stdTrsGeoBuilder;
  
  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> m_convergenceMtd;

  /// variable set in which linearizing the residual
  Common::SelfRegistPtr<Framework::ConvectiveVarSet>  m_linearVar;

  /// distribution data
  DistributionData m_distData;
  
  /// residual factor
  CFreal m_resFactor;
  
  /// are we in the middle of the initialization phase
  bool m_isInitializationPhase;

  /// string for configuration of the scalar splitter
  std::string m_scalarSplitterStr;

  /// string for configuration of the system splitter for implict computation; used to compute the Jacobian
  std::string m_jacobSysSplitterStr;

  /// string for configuration of the scalar splitter for implict computation; used to compute the Jacobian
  std::string m_jacobSclSplitterStr;

  /// string for configuration of the source term splitter
  std::vector<std::string> m_sourceTermSplitterStr;

  /// string for configuration of the distribution variables
  std::string m_distribVarStr;

  /// string for configuration of the linearizing variables
  std::string m_linearVarStr;

  /// string for configuration of the source term computer
  std::vector<std::string> m_stComputerStr;

  /// string for configuration of the diffusive term computer
  std::string m_diffTermComputerStr;

  /// string for configuration of the jacobian fix computer
  std::string m_jacobFixComputerStr;

  /// string for configuration of the linearizer
  std::string m_linearizerStr;

  /// string for configuration of the transformation from solution to
  /// distribution vars
  std::string m_solutionToDistMatTransStr;

  /// string for configuration of the transformation from distribution
  /// to solution vars
  std::string m_distToSolutionMatTransStr;

  /// string for configuration of the transformation from linear
  /// to distribution vars
  std::string m_linearToDistMatTransStr;

  /// string for configuration of the transformation from linear
  /// to solution in solution vars
  std::string m_solutionToLinearInUpdateMatTransStr;

  /// string for configuration of the transformation from solution
  /// to linear vars
  std::string m_solutionToLinearMatTransStr;

  /// string for configuration of the transformation from update
  /// to linear vars
  std::string m_updateToLinearVecTransStr;

  /// string for configuration of the transformation from update
  /// to solution vars
  std::string m_updateToSolutionVecTransStr;

  /// string for configuration of the transformation from solution
  /// to update in update vars
  std::string m_solToUpdateInUpdateMatTransStr;

  /// string for configuration of the transformation from update to
  /// solution in update vars
  std::string m_updateToSolutionInUpdateMatTransStr;

  /// string for the configuration of the numerical integrator Order
  std::string m_integratorOrderStr;

  /// string for the configuration of the numerical contour integrator Order
  std::string m_cIntegratorOrderStr;

  /// string for the configuration of the numerical volume integrator Order
  std::string m_vIntegratorOrderStr;

  /// string for the configuration of the numerical integrator QuadratureType
  std::string m_integratorQuadratureStr;

  /// flag that tells if the SubSystem is axisymmetric
  bool m_isAxisymm;

  /// flag that tells if a source term has to be included in the convective flux
  bool m_includeSourceInFlux;

  /// Fluctuation splitting strategy name
  std::string m_fsStrategyName;
  /// Fluctuation splitting strategy name
  std::string m_adStrategyName;
  /// Fluctuation Split Strategy
  Common::SelfRegistPtr<FluctuationSplitStrategy> m_fsStrategy;

  /// Fluctuation Split Strategy
  Common::SelfRegistPtr<ArtificialDiffusionStrategy> m_adStrategy;

  /// Jacobian computation strategy name
  std::string m_jacobStrategyName;
  
  /// Fluctuation Split Strategy
  Common::SelfRegistPtr<ComputeJacobStrategy> m_jacobStrategy;

  /// The contour Integrator
  Framework::ContourIntegrator m_ContourIntegrator;

  /// The VolumeIntegrator
  Framework::VolumeIntegrator m_VolumeIntegrator;

  /// add an artificial diffusion term?
  bool m_artdiff;
  
  /// scalar first
  bool m_scalarFirst;

}; // end of class FluctuationSplitData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for FluctSplit
typedef Framework::MethodCommand<FluctuationSplitData> FluctuationSplitCom;

/// Definition of a command provider for FluctSplit
typedef Framework::MethodCommand<FluctuationSplitData>::PROVIDER FluctuationSplitComProvider;

/// Definition of a command for FluctSplit
typedef Framework::MethodStrategy<FluctuationSplitData> FluctuationSplitStrat;

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_FluctuationSplitData_hh
