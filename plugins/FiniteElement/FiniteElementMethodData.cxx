


#include "Common/NullPointerException.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FiniteElement/ComputeConvectiveTerm.hh"
#include "FiniteElement/ComputeDiffusiveTerm.hh"
#include "FiniteElement/ComputeLinearSourceTerm.hh"
#include "FiniteElement/ComputeIndepSourceTerm.hh"
#include "FiniteElement/ComputeInertiaTerm.hh"
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/FiniteElementMethodData.hh"
#include "FiniteElement/ComputeResidualStrategy.hh"
#include "FiniteElement/ComputeJacobStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<FiniteElementMethodData>, FiniteElementMethodData, FiniteElementModule> nullFiniteElementMethodComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethodData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("IntegratorOrder","Order of the Integration to be used for numerical quadrature.");
   options.addConfigOption< std::string >("JacobianStrategy","Strategy to compute the jacobian system matrix.");
   options.addConfigOption< std::string >("SourceVar","Source variable set.");
   options.addConfigOption< std::string >("IndepSourceEntity","Source entity.");
   options.addConfigOption< std::string >("InertiaVar","Inertia variable set.");
   options.addConfigOption< std::string >("ResidualStrategy","Strategy to compute the system residual.");
   options.addConfigOption< std::string >("IntegratorQuadrature","Type of Quadrature to be used in the Integration.");
}

//////////////////////////////////////////////////////////////////////////////

FiniteElementMethodData::FiniteElementMethodData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  _convergenceMtd(),
  _inertiaVar(),
  _sourceVar(),
  _convectiveTerm(),
  _diffusiveTerm(),
  _inertiaTerm(),
  _linearSourceTerm(),
  _indepSourceTerm(),
  _jacobianStrategy(),
  _residualStrategy(),
  _stdTrsGeoBuilder(),
  _resFactor(1.0),
  m_local_elem_data()
{
   addConfigOptionsTo(this);

  _inertiaVarStr = "Null";
   setParameter("InertiaVar",&_inertiaVarStr);

  _sourceVarStr = "Null";
   setParameter("SourceVar",&_sourceVarStr);

  _indepSourceEntityStr = "Galerkin";
   setParameter("IndepSourceEntity",&_indepSourceEntityStr);

  _integratorOrderStr = "P1";
   setParameter("IntegratorOrder",&_integratorOrderStr);

  _integratorQuadratureStr = "INVALID";
   setParameter("IntegratorQuadrature",&_integratorQuadratureStr);

   _jacobianStrategyStr = "Numerical";
   setParameter("JacobianStrategy",&_jacobianStrategyStr);

   _residualStrategyStr = "StdElementComputer";
   setParameter("ResidualStrategy",&_residualStrategyStr);
}

//////////////////////////////////////////////////////////////////////////////

FiniteElementMethodData::~FiniteElementMethodData()
{
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethodData::configure ( Config::ConfigArgs& args )
{
  SpaceMethodData::configure(args);

  /// @todo add here the setup for the specific integrators for each element
  ///       that has a different set of shape and interpolator type.

  configureVarSets(args);

  configureTerms(args);

  configureIntegrator(args);

  SharedPtr<FiniteElementMethodData> thisPtr(this);

    CFLogDebugMax( "Configure Strategy Type: " << _residualStrategyStr << "\n");

    try {
      Common::SafePtr<
      BaseMethodStrategyProvider<FiniteElementMethodData,ComputeResidualStrategy> > prov =
        Environment::Factory<ComputeResidualStrategy>::getInstance().getProvider(_residualStrategyStr);
      cf_assert(prov.isNotNull());
      _residualStrategy = prov->create(_residualStrategyStr,thisPtr);
      configureNested ( _residualStrategy.getPtr(), args );
    }
    catch (Common::NoSuchValueException& e) {
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing Null of Type: "
                  << ComputeResidualStrategy::getClassName()
                  << " instead ...\n");
      Common::SafePtr<
      BaseMethodStrategyProvider<FiniteElementMethodData,ComputeResidualStrategy> > prov =
        Environment::Factory<ComputeResidualStrategy>::getInstance().getProvider("Null");
      cf_assert(prov.isNotNull());
      _residualStrategy = prov->create("Null", thisPtr);
    }
    cf_assert(_residualStrategy.isNotNull());

    CFLogDebugMax( "Configure Strategy Type: " << _jacobianStrategyStr << "\n");

    try {
      Common::SafePtr<
      BaseMethodStrategyProvider<FiniteElementMethodData,ComputeJacobStrategy> > prov =
        Environment::Factory<ComputeJacobStrategy>::getInstance().getProvider(_jacobianStrategyStr);
      cf_assert(prov.isNotNull());
      _jacobianStrategy = prov->create(_jacobianStrategyStr,thisPtr);
      configureNested ( _jacobianStrategy.getPtr(), args );
    }
    catch (Common::NoSuchValueException& e) {
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing Null of Type: "
                  << ComputeJacobStrategy::getClassName()
                  << " instead ...\n");
      Common::SafePtr<
      BaseMethodStrategyProvider<FiniteElementMethodData,ComputeJacobStrategy> > prov =
        Environment::Factory<ComputeJacobStrategy>::getInstance().getProvider("Null");
      cf_assert(prov.isNotNull());
      _jacobianStrategy = prov->create("Null",thisPtr);
    }
    cf_assert(_jacobianStrategy.isNotNull());

  CFLog(INFO,"FiniteElementMethod: Using Integrator Quadrature: "
             << _integratorQuadratureStr << "\n");
  CFLog(INFO,"FiniteElementMethod: Using Integrator Order: "
             << _integratorOrderStr << "\n");
  CFLog(INFO,"FiniteElementMethod: Using Convective VarSet: "
             << _updateVarStr << "\n");
  CFLog(INFO,"FiniteElementMethod: Using Diffusive VarSet: "
             << _diffusiveVarStr << "\n");
  CFLog(INFO,"FiniteElementMethod: Using Inertia VarSet: "
             << _inertiaVarStr << "\n");
  CFLog(INFO,"FiniteElementMethod: Using Source VarSet: "
             << _sourceVarStr << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethodData::configureTerms ( Config::ConfigArgs& args )
{

  const std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = 
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  //by default we dont have convection
  ///@todo change this when we will start using ConvTerms!!
  std::string convectiveTermStr = "NullConvectiveTerm";
  std::string convectiveEntityStr = "Galerkin";

  std::string diffusiveTermStr = "NullDiffusiveTerm";
  std::string diffusiveEntityStr = "Null";
  if(getDiffusiveVar()->getName() != "Null"){
    diffusiveTermStr = "DiffusiveTerm";
    diffusiveEntityStr = "Galerkin" +
                         getDiffusiveVar()->getName();
  }

  //by default we dont have Inertia
  std::string inertiaTermStr = "NullInertiaTerm";
  std::string inertiaEntityStr = "Null";
  if(_inertiaVarStr != "Null"){
    inertiaTermStr = "InertiaTerm";
    inertiaEntityStr = "Galerkin"  +
                       getInertiaVar()->getName();;
  }

  //by default we dont have Linear Source Term
  std::string linearSourceTermStr = "NullLinearSourceTerm";
  std::string linearSourceEntityStr = "Galerkin";
  if(_sourceVarStr != "Null") linearSourceTermStr = "LinearSourceTerm";

  //by default we dont have Indep Source Term
  std::string indepSourceTermStr = "NullIndepSourceTerm";
  if(_sourceVarStr != "Null") indepSourceTermStr = "IndepSourceTerm";

  SharedPtr<FiniteElementMethodData> thisPtr(this);

//   _convectiveTerm.reset(new ComputeConvectiveTerm("ConvectiveTerm"));
//   _diffusiveTerm.reset(new ComputeDiffusiveTerm("DiffusiveTerm"));
//   _inertiaTerm.reset(new ComputeInertiaTerm("InertiaTerm"));
//   _linearSourceTerm.reset(new ComputeLinearSourceTerm("LinearSourceTerm"));
//   _indepSourceTerm.reset(new ComputeIndepSourceTerm("IndepSourceTerm"));

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,ComputeConvectiveTerm> > provConv =
    Environment::Factory<ComputeConvectiveTerm>::getInstance().getProvider(convectiveTermStr);
  cf_assert(provConv.isNotNull());
  _convectiveTerm = provConv->create(convectiveTermStr,thisPtr);
  configureNested ( _convectiveTerm.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,ConvectiveEntity> > provConvEnt =
    Environment::Factory<ConvectiveEntity>::getInstance().getProvider(convectiveEntityStr);
  cf_assert(provConvEnt.isNotNull());
  _convectiveEntity = provConvEnt->create(convectiveEntityStr,thisPtr);
  configureNested ( _convectiveEntity.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,ComputeDiffusiveTerm> > provDiff =
    Environment::Factory<ComputeDiffusiveTerm>::getInstance().getProvider(diffusiveTermStr);
  cf_assert(provDiff.isNotNull());
  _diffusiveTerm = provDiff->create(diffusiveTermStr,thisPtr);
  configureNested ( _diffusiveTerm.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,DiffusiveEntity> > provDiffEnt =
    Environment::Factory<DiffusiveEntity>::getInstance().getProvider(diffusiveEntityStr);
  cf_assert(provDiffEnt.isNotNull());
  _diffusiveEntity = provDiffEnt->create(diffusiveEntityStr,thisPtr);
  configureNested ( _diffusiveEntity.getPtr(), args );


  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,ComputeInertiaTerm> > provInertia =
    Environment::Factory<ComputeInertiaTerm>::getInstance().getProvider(inertiaTermStr);
  cf_assert(provInertia.isNotNull());
  _inertiaTerm = provInertia->create(inertiaTermStr,thisPtr);
  configureNested ( _inertiaTerm.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,InertiaEntity> > provInertiaEnt =
    Environment::Factory<InertiaEntity>::getInstance().getProvider(inertiaEntityStr);
  cf_assert(provInertiaEnt.isNotNull());
  _inertiaEntity = provInertiaEnt->create(inertiaEntityStr,thisPtr);
  configureNested ( _inertiaEntity.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,ComputeLinearSourceTerm> > provLinear =
    Environment::Factory<ComputeLinearSourceTerm>::getInstance().getProvider(linearSourceTermStr);
  cf_assert(provLinear.isNotNull());
  _linearSourceTerm = provLinear->create(linearSourceTermStr,thisPtr);
  configureNested ( _linearSourceTerm.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,LinearSourceEntity> > provLinearSourceEnt =
    Environment::Factory<LinearSourceEntity>::getInstance().getProvider(linearSourceEntityStr);
  cf_assert(provLinearSourceEnt.isNotNull());
  _linearSourceEntity = provLinearSourceEnt->create(linearSourceEntityStr,thisPtr);
  configureNested ( _linearSourceEntity.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,ComputeIndepSourceTerm> > provIndep =
    Environment::Factory<ComputeIndepSourceTerm>::getInstance().getProvider(indepSourceTermStr);
  cf_assert(provIndep.isNotNull());
  _indepSourceTerm = provIndep->create(indepSourceTermStr,thisPtr);
  configureNested ( _indepSourceTerm.getPtr(), args );

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,IndepSourceEntity> > provIndepSourceEnt =
    Environment::Factory<IndepSourceEntity>::getInstance().getProvider(_indepSourceEntityStr);
  cf_assert(provIndepSourceEnt.isNotNull());
  _indepSourceEntity = provIndepSourceEnt->create(_indepSourceEntityStr,thisPtr);
  configureNested ( _indepSourceEntity.getPtr(), args );

}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethodData::configureIntegrator ( Config::ConfigArgs& args )
{
  const CFQuadrature::Type quadType =
    CFQuadrature::Convert::to_enum( _integratorQuadratureStr );

  const CFPolyOrder::Type order =
    CFPolyOrder::Convert::to_enum( _integratorOrderStr );

  _femVolumeIntegrator.setIntegrationForAllGeo(quadType,order);

}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethodData::configureVarSets ( Config::ConfigArgs& args )
{
  /// @todo Framework::VarSet is not a ConfigObject. should it be?

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
//   // Update VarSet
//   CFLog(INFO,"Configuring Convective VarSet: " <<_updateVarStr<<" \n");
//   _updateVar.reset(Environment::Factory<VarSet>::getInstance().getProvider(_updateVarStr)->
//     create(physModel->getImplementor()->getConvectiveTerm()) );
//   cf_assert(_updateVar.isNotNull());
// //  configureNested ( _updateVar.getPtr(), args );
//
//   // Diffusive VarSet
//   CFLog(INFO,"Configuring Diffusive VarSet: " <<_diffusiveVarStr<<" \n");
//   _diffusiveVar.reset(Environment::Factory<DiffusiveVarSet>::getInstance().getProvider(_diffusiveVarStr)->
//     create(_diffusiveVarStr, physModel->getImplementor()) );
//   cf_assert(_diffusiveVar.isNotNull());
//   configureNested ( _diffusiveVar.getPtr(), args );

  CFLog(INFO,"Configuring Inertia VarSet: " <<_inertiaVarStr<<" \n");
  // Inertia VarSet
  _inertiaVar.reset(Environment::Factory<InertiaVarSet>::getInstance().getProvider(_inertiaVarStr)->
    create(_inertiaVarStr) );
  cf_assert(_inertiaVar.isNotNull());
  configureNested ( _inertiaVar.getPtr(), args );

  // Source VarSet
  CFLog(INFO,"Configuring Source VarSet: " <<_sourceVarStr<<" \n");
  _sourceVar.reset(Environment::Factory<SourceVarSet>::getInstance().getProvider(_sourceVarStr)->
    create(_sourceVarStr));
  cf_assert(_sourceVar.isNotNull());
  configureNested ( _sourceVar.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethodData::setup()
{
  CFAUTOTRACE;

  SpaceMethodData::setup();

  // set up the volume integrator
  getFEMVolumeIntegrator()->setup();

//   CFLog(VERBOSE,"Setting up VarSets \n");

  _updateVar->setup();
  _diffusiveVar->setup();
  _inertiaVar->setup();
  _sourceVar->setup();

//   CFLog(VERBOSE,"Setting up Terms\n");
  _convectiveTerm->setup();
  _diffusiveTerm->setup();
  _inertiaTerm->setup();
  _linearSourceTerm->setup();
  _indepSourceTerm->setup();

  _diffusiveEntity->setup();
  _convectiveEntity->setup();
  _inertiaEntity->setup();
  _linearSourceEntity->setup();
  _indepSourceEntity->setup();

  // set up the GeometricEntity builder
  _stdTrsGeoBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

void FiniteElementMethodData::unsetup()
{
  _stdTrsGeoBuilder.unsetup();

  getFEMVolumeIntegrator()->unsetup();

  SpaceMethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<VolumeIntegrator> FiniteElementMethodData::getVolumeIntegrator()
{
  return &_femVolumeIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<FEM_VolumeIntegrator> FiniteElementMethodData::getFEMVolumeIntegrator()
{
  return &_femVolumeIntegrator;
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

