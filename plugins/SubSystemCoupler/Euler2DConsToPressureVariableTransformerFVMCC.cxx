#include "Euler2DConsToPressureVariableTransformerFVMCC.hh"
#include "Framework/PhysicalModel.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "SubSysCouplerData.hh"

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

MethodStrategyProvider<Euler2DConsToPressureVariableTransformerFVMCC,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
euler2DConsToPressureVariableTransformerFVMCCProvider("Euler2DConsToPressureFVMCC");

/////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPressureVariableTransformerFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("ReferencePressure","Value of pressure to substract from pressure of the flow.");
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPressureVariableTransformerFVMCC::Euler2DConsToPressureVariableTransformerFVMCC(const std::string& name) :
  Euler2DConsToPressureVariableTransformer(name)
{
   addConfigOptionsTo(this);

   _referencePressure = 0.;
   setParameter("ReferencePressure",&_referencePressure);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPressureVariableTransformerFVMCC::~Euler2DConsToPressureVariableTransformerFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector* Euler2DConsToPressureVariableTransformerFVMCC::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original)
{
  cf_assert(original.size() == 4);

  const CFreal rho  = original[0];
  const CFreal rhoU = original[1];
  const CFreal rhoV = original[2];
  const CFreal rhoE = original[3];
  const CFreal V2 = (rhoU*rhoU + rhoV*rhoV)/(rho*rho);
  const CFreal localPressure = (_model->getConvTerm()->getGamma() - 1.)*(rhoE - 0.5*rho*V2);
  const CFreal pressure = localPressure - _referencePressure;

  /// builder for standard TRS GeometricEntity's
  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
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

  const CFreal norm = normal.norm2();
  normal[XX] /= norm;
  normal[YY] /= norm;

  _transVector[0] = pressure * normal[XX];
  _transVector[1] = pressure * normal[YY];

  //release GeometricEntity
  geoBuilder->releaseGE();

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
