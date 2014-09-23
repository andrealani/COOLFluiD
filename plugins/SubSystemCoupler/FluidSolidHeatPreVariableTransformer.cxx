#include "FluidSolidHeatPreVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethodData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FluidSolidHeatPreVariableTransformer,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
FluidSolidHeatPreVariableTransformerProvider("FluidSolidHeatPre");

/////////////////////////////////////////////////////////////////////////////

void FluidSolidHeatPreVariableTransformer::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("hValue","Value given to h.");
}

//////////////////////////////////////////////////////////////////////

FluidSolidHeatPreVariableTransformer::FluidSolidHeatPreVariableTransformer(const std::string& name) :
  PreVariableTransformer(name),
  socket_faceNeighCell("faceNeighCell"),
  _diffVar(CFNULL)
{
   addConfigOptionsTo(this);

   _hConst = 3000.;
   setParameter("hValue",&_hConst);

}

//////////////////////////////////////////////////////////////////////

FluidSolidHeatPreVariableTransformer::~FluidSolidHeatPreVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FluidSolidHeatPreVariableTransformer::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = PreVariableTransformer::needsSockets();

  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////

void FluidSolidHeatPreVariableTransformer::setup()
{

  PreVariableTransformer::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  _transVector.resize(getTransformedSize(2));
  _gradients.resize(nbEqs);

  for(CFuint iEq=0; iEq < nbEqs; iEq++)
  {
    _gradients[iEq]->resize(nbDim);
  }
  _coords.resize(1);
  _coords[0].resize(nbDim);
  _normal.resize(nbDim);

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  cf_assert(spaceMethod.isNotNull());

  _diffVar = spaceMethod->getSpaceMethodData()->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  cf_assert(_diffVar.isNotNull());

}

//////////////////////////////////////////////////////////////////////

RealVector* FluidSolidHeatPreVariableTransformer::preTransform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original)
{
  cf_assert(original.size() == 4);
  cf_assert(_coordType == "Gauss");
  cf_assert(faces.size() == 1);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  /// builder for standard TRS GeometricEntity's
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  ///Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  GeometricEntity* currFace = geoBuilder->buildGE();
  const CFuint faceID = currFace->getID();

  ///release GeometricEntity
  geoBuilder->releaseGE();

  /// Build the neighbor cell
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;

  geoData.trs = cellTrs;
  geoData.idx = cellTrsID;
  GeometricEntity* neighborCell = geoBuilder->buildGE();

  ///@todo modify computeShapeFunctionGradients to
  ///accept RealVector and not only vector<RealVector>
  _coords[0] = coord;

  /// Compute the shape function gradient at the projected point
  std::vector<RealMatrix> cellGradients = neighborCell->computeSolutionShapeFunctionGradients(_coords);

  const CFuint nbNodesInCell = neighborCell->getNodes()->size();

  /// Compute the gradients of T
  for (CFuint iEq = 0; iEq < nbEqs ; ++iEq){
    *(_gradients[iEq]) = 0.;
  }

  for (CFuint iEq = 0; iEq < nbEqs ; ++iEq){
    for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode){
      for (CFuint iDim =0; iDim < nbDim ; ++iDim){
        (*(_gradients[iEq]))[iDim] += (cellGradients[0])(iNode,iDim) * (*(neighborCell->getState(iNode)))[iEq] ;
      }
    }
  }

  /// Get the face normal
  ///@todo change this to use not the
  /// Average but the local FaceNormal (for 2nd order elements)
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  _normal = (neighborCell->computeAvgFaceNormals())[iFaceLocal];
  _normal.normalize();

  /// Project the gradient of T on the face normal
  ///@here we assume that we are in Puvt
  CFreal dTdn = 0.;
  for (CFuint iDim=0; iDim < nbDim ; ++iDim){
    dTdn += (*(_gradients[3]))[iDim] * _normal[iDim];
  }

  ///release GeometricEntity
  geoBuilder->releaseGE();

  /// flux = lambda * gradientT

  //First compute the thermal conductivity
  ///This is a transformer for heat transfer between fluid and solid
  ///We are at the interface, so at the wall => wall distance is 0.
  _diffVar->setWallDistance(0.);
  const CFreal mu = _diffVar->getDynViscosity(original, _gradients);
  const CFreal lambda = _diffVar->getThermConductivity(original, mu);

  ///This is what T. Verstraeten does in his article asme GT2006-90161
  const CFreal currentFlux = lambda*dTdn;
  ///@todo here we suppose that we use puvT variables check that it is ok
  const CFreal currentTemperature = original[3];
  const CFreal Tfl = currentTemperature - currentFlux/_hConst;

  _transVector[0] = _hConst;
  _transVector[1] = Tfl;

  return (&_transVector);

}

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
