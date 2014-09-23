#include "Euler2DConsToPressurePredictorVariableTransformerFVMCC.hh"
#include "Framework/PhysicalModel.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "SubSysCouplerData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Euler2DConsToPressurePredictorVariableTransformerFVMCC,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
Euler2DConsToPressurePredictorVariableTransformerFVMCCProvider("Euler2DConsToPressurePredictorFVMCC");

/////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPressurePredictorVariableTransformerFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("ReferencePressure","Value of pressure to substract from pressure of the flow.");
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPressurePredictorVariableTransformerFVMCC::Euler2DConsToPressurePredictorVariableTransformerFVMCC(const std::string& name) :
  Euler2DConsToPressureVariableTransformer(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  _pastPressures(0)
{
   addConfigOptionsTo(this);

   _referencePressure = 0.;
   setParameter("ReferencePressure",&_referencePressure);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPressurePredictorVariableTransformerFVMCC::~Euler2DConsToPressurePredictorVariableTransformerFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPressurePredictorVariableTransformerFVMCC::setup()
{
  PreVariableTransformer::setup();

  _transVector.resize(getTransformedSize(4));
}


//////////////////////////////////////////////////////////////////////////////

RealVector* Euler2DConsToPressurePredictorVariableTransformerFVMCC::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original)
{
  cf_assert(original.size() == 4);

  const CFreal rho  = original[0];
  const CFreal rhoU = original[1];
  const CFreal rhoV = original[2];
  const CFreal rhoE = original[3];
  const CFreal V2 = (rhoU*rhoU + rhoV*rhoV)/(rho*rho);

  const CFreal localPressure = (_model->getConvTerm()->getGamma() - 1.)*(rhoE - 0.5*rho*V2);
  const CFreal pressure = localPressure - _referencePressure;

  if(SubSystemStatusStack::getActive()->isSubIterationFirstStep() &&
      (SubSystemStatusStack::getActive()->getNbIter() <= 1))
  {
    if(_iOtherState >= _pastPressures.size()){
      _pastPressures.resize(_iOtherState+1);
    }
  }

  if(SubSystemStatusStack::getActive()->getNbIter() == 1)
  {
    _pastPressures[_iOtherState] = pressure;
  }

  const CFreal pastPressure = _pastPressures[_iOtherState];
  const CFreal predictedPressure = (2. * pressure) - pastPressure;

  if((SubSystemStatusStack::getActive()->doingSubIterations() &&
     SubSystemStatusStack::getActive()->isSubIterationLastStep()) ||
     (!SubSystemStatusStack::getActive()->doingSubIterations()))
  {
    cf_assert(_iOtherState < _pastPressures.size());
    _pastPressures[_iOtherState] = pressure;
  }

  /// builder for standard TRS GeometricEntity's
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  geoData.isBFace = true;

  GeometricEntity* currFace = geoBuilder->buildGE();

  /// Get the face normal
  ///@todo change this to use not the
  /// Average but the local FaceNormal (for 2nd order elements)
  RealVector normal = currFace->computeAvgCellNormal();

  CFreal norm = normal.norm2();
  normal[XX] /= norm;
  normal[YY] /= norm;

  _transVector[0] = predictedPressure * normal[XX];
  _transVector[1] = predictedPressure * normal[YY];

  ///release GeometricEntity
  geoBuilder->releaseGE();

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
