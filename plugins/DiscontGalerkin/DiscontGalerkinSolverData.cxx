

#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"
//#include "DiscontGalerkin/DiscontGalerkinStrategy.hh"

#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"
// #include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<
  NullMethodCommand< DiscontGalerkinSolverData >,DiscontGalerkinSolverData,DiscontGalerkinModule >
  nullDiscontGalerkinSolverComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolverData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("VolumeIntegratorOrder","Order of the Integration to be used for numerical quadrature for comutation of volume integrals.");
  options.addConfigOption< std::string >("VolumeIntegratorQuadrature","Type of Quadrature to be used in the Integration for comutation of volume integrals.");
  options.addConfigOption< std::string >("ContourIntegratorOrder","Order of the Integration to be used for numerical quadrature for comutation of contour integrals.");
  options.addConfigOption< std::string >("ContourIntegratorQuadrature","Type of Quadrature to be used in the Integration for comutation of contour integrals.");
  options.addConfigOption< std::string >("Sigma","Sigma constant in computation of viscous flow");
  options.addConfigOption< std::string >("TypeOfDGFEM","Theta constant in computation of viscous flow.");
  options.addConfigOption< std::string >("Reynolds","Reynolds constant in computation of viscous flow.");
  options.addConfigOption< std::string >("MaxCFL","Constant in computation of time step.");
  options.addConfigOption< std::string >("Alpha","Constant in computation of time step.");
//   options.addConfigOption< std::string >("StrategyForSomething","A MethodStrategy to be used for some calculation (default = DiscontGalerkinStrategy).");
}

//////////////////////////////////////////////////////////////////////////////

DiscontGalerkinSolverData::DiscontGalerkinSolverData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder(),
  m_FaceBuilder(),
  m_maxEigenValue(0),
  m_Residual(0)
{
  addConfigOptionsTo(this);

  m_volintorderStr = "P1";
  setParameter( "VolumeIntegratorOrder",      &m_volintorderStr );

  m_volintquadStr  = "INVALID";
  setParameter( "VolumeIntegratorQuadrature", &m_volintquadStr );

  m_conintorderStr = "P1";
  setParameter( "ContourIntegratorOrder",      &m_conintorderStr );

  m_conintquadStr  = "INVALID";
  setParameter( "ContourIntegratorQuadrature", &m_conintquadStr );

  m_sigmaStr  = "100";
  setParameter( "Sigma", &m_sigmaStr );

  m_typeStr  = "";
  setParameter( "TypeOfDGFEM", &m_typeStr );

  m_reynoldStr  = "";
  setParameter( "Reynolds", &m_reynoldStr );

  m_maxCFLStr  = "100 000";
  setParameter( "MaxCFL", &m_maxCFLStr );

  m_alphaStr  = "1.0";
  setParameter( "Alpha", &m_alphaStr );

//   m_emptystrategyStr = "DiscontGalerkinStrategy";
//   setParameter( "StrategyForSomething", &m_emptystrategyStr );
}

//////////////////////////////////////////////////////////////////////////////

DiscontGalerkinSolverData::~DiscontGalerkinSolverData()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolverData::configure ( Config::ConfigArgs& args )
{
  SpaceMethodData::configure(args);
  SharedPtr< DiscontGalerkinSolverData > thisPtr(this);

  configureIntegrator();

  /* add here the setup for the specific integrators for each element that
   * have different set of shape and interpolator type */

  CFLog(INFO,"DiscontGalerkinSolver: volume integrator quadrature: " << m_volintquadStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: volume integrator order: " << m_volintorderStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: contour integrator quadrature: " << m_conintquadStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: contour integrator order: " << m_conintorderStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: constant sigma: " << m_sigmaStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: constant type of DGFEM: " << m_typeStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: constant reynolds: " << m_reynoldStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: constant maxCFL: " << m_maxCFLStr << "\n");
  CFLog(INFO,"DiscontGalerkinSolver: constant alpha: " << m_alphaStr << "\n");
  /* add here different strategies configuration */

/*  CFLog(INFO,"Configure strategy type: " << m_emptystrategyStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< DiscontGalerkinSolverData, DiscontGalerkinStrategy > >
      prov = Environment::Factory< DiscontGalerkinStrategy >::getInstance().getProvider(
        m_emptystrategyStr );
    cf_assert(prov.isNotNull());
    m_emptystrategy = prov->create(m_emptystrategyStr,thisPtr);
    configureNested ( m_emptystrategy.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  DiscontGalerkinStrategy ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< DiscontGalerkinSolverData, DiscontGalerkinStrategy > >
      prov = Environment::Factory< DiscontGalerkinStrategy >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_emptystrategy = prov->create("Null", thisPtr);

  }
  cf_assert(m_emptystrategy.isNotNull());
*/
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolverData::setup()
{
  CFAUTOTRACE;

  // setup TRS Geo builder
  m_stdTrsGeoBuilder.setup();

  // setup face builder
  m_FaceBuilder.setup();

  m_volumeIntegrator.setup();
  m_contourIntegrator.setup();
}

//////////////////////////////////////////////////////////////////////////////

void DiscontGalerkinSolverData::configureIntegrator()
{

  CFQuadrature::Type quadType = CFQuadrature::Convert::to_enum (m_volintquadStr);

  CFPolyOrder::Type order = CFPolyOrder::Convert::to_enum (m_volintorderStr);

  m_volumeIntegrator.setIntegrationForAllGeo(quadType,order);


  quadType =  CFQuadrature::Convert::to_enum (m_conintquadStr);

  order = CFPolyOrder::Convert::to_enum (m_conintorderStr);

  m_contourIntegrator.setIntegrationForAllGeo(quadType,order);
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<VolumeIntegrator> DiscontGalerkinSolverData::getVolumeIntegrator()
{
  return &m_volumeIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<ContourIntegrator> DiscontGalerkinSolverData::getContourIntegrator()
{
  return &m_contourIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

