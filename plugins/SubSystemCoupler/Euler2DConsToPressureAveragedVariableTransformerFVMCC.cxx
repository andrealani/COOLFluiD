#include "Euler2DConsToPressureAveragedVariableTransformerFVMCC.hh"
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

MethodStrategyProvider<Euler2DConsToPressureAveragedVariableTransformerFVMCC,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
euler2DConsToPressureAveragedVariableTransformerFVMCCProvider("Euler2DConsToPressureAveragedFVMCC");

/////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPressureAveragedVariableTransformerFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("ReferencePressure","Value of pressure to substract from pressure of the flow.");
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPressureAveragedVariableTransformerFVMCC::Euler2DConsToPressureAveragedVariableTransformerFVMCC(const std::string& name) :
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

Euler2DConsToPressureAveragedVariableTransformerFVMCC::~Euler2DConsToPressureAveragedVariableTransformerFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPressureAveragedVariableTransformerFVMCC::setup()
{
  PreVariableTransformer::setup();

  _transVector.resize(getTransformedSize(4));
}

//////////////////////////////////////////////////////////////////////////////

RealVector* Euler2DConsToPressureAveragedVariableTransformerFVMCC::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original)
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
      (SubSystemStatusStack::getActive()->getNbIter() == 1))
  {
    if(_iOtherState+1 > _pastPressures.size()){
      _pastPressures.resize(_iOtherState+1);
    }
  }

  if(SubSystemStatusStack::getActive()->getNbIter() == 1)
  {
    _pastPressures[_iOtherState] = pressure;
  }

  const CFreal pastPressure = _pastPressures[_iOtherState];
  const CFreal avgPressure = (pressure + pastPressure)*0.5;

  if(SubSystemStatusStack::getActive()->doingSubIterations() &&
     SubSystemStatusStack::getActive()->isSubIterationLastStep())
  {
    cf_assert(_iOtherState < _pastPressures.size());
    _pastPressures[_iOtherState] = pressure;
  }

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> geoBuilder;
  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder.getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  geoBuilder.setup();

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  geoData.isBFace = true;

  GeometricEntity& currFace = *geoBuilder.buildGE();

  /// Get the face normal
  ///@todo change this to use not the
  /// Average but the local FaceNormal (for 2nd order elements)
  RealVector normal = currFace.computeAvgCellNormal();

  CFreal norm = normal.norm2();
  normal[XX] /= norm;
  normal[YY] /= norm;

  _transVector[0] = avgPressure * normal[XX];
  _transVector[1] = avgPressure * normal[YY];

  ///release GeometricEntity
  geoBuilder.releaseGE();

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
