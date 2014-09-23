#include "Euler2DConsToPressureVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Euler2DConsToPressureVariableTransformer,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
euler2DConsToPressureVariableTransformerProvider("Euler2DConsToPressure");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPressureVariableTransformer::Euler2DConsToPressureVariableTransformer(const std::string& name) :
  PreVariableTransformer(name),
  _model(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPressureVariableTransformer::~Euler2DConsToPressureVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPressureVariableTransformer::setup()
{

  PreVariableTransformer::setup();

  _model = Framework::PhysicalModelStack::getActive()->getImplementor().d_castTo<EulerPhysicalModel<DIM_2D> >();
  _transVector.resize(getTransformedSize(4));
}


//////////////////////////////////////////////////////////////////////////////

RealVector* Euler2DConsToPressureVariableTransformer::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original)
{
  cf_assert(original.size() == 4);

  const CFreal rho  = original[0];
  const CFreal rhoU = original[1];
  const CFreal rhoV = original[2];
  const CFreal rhoE = original[3];
  const CFreal V2 = (rhoU*rhoU + rhoV*rhoV)/(rho*rho);

  const CFreal pressure = (_model->getConvTerm()->getGamma() - 1.)*(rhoE - 0.5*rho*V2);

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  GeometricEntity& currFace = *geoBuilder.buildGE();

  /// Get the face normal
  ///@todo change this to use not the
  /// Average but the local FaceNormal (for 2nd order elements)
  RealVector normal = currFace.computeAvgCellNormal();

  CFreal norm = normal.norm2();
  normal[XX] /= norm;
  normal[YY] /= norm;

  _transVector[0] = pressure * normal[XX];
  _transVector[1] = pressure * normal[YY];

  ///release GeometricEntity
  geoBuilder.releaseGE();

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
