#include "Temp2TempAndFluxVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "SubSystemCoupler/SubSystemCouplerHeat.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Heat/HeatPhysicalModel.hh"
#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::Heat;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Temp2TempAndFluxVariableTransformer,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerHeatModule>
temp2TempAndFluxVariableTransformerProvider("Temp2TempAndFlux");

//////////////////////////////////////////////////////////////////////////////

Temp2TempAndFluxVariableTransformer::Temp2TempAndFluxVariableTransformer(const std::string& name) :
  PreVariableTransformer(name),
  _model(CFNULL),
  socket_faceNeighCell("faceNeighCell")
{
}

//////////////////////////////////////////////////////////////////////////////

Temp2TempAndFluxVariableTransformer::~Temp2TempAndFluxVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Temp2TempAndFluxVariableTransformer::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Temp2TempAndFluxVariableTransformer::setup()
{

  PreVariableTransformer::setup();

  _model = Framework::PhysicalModelStack::getActive()->getImplementor().d_castTo<HeatPhysicalModel>();
  _transVector.resize(getTransformedSize(PhysicalModelStack::getActive()->getNbEq()));
  _gradientsT.resize(PhysicalModelStack::getActive()->getDim());
  _coords.resize(1);
  _coords[0].resize(PhysicalModelStack::getActive()->getDim());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

RealVector* Temp2TempAndFluxVariableTransformer::preTransform(const vector<GeoEntityIdx>& faces, const RealVector& coord,const RealVector& original)
{
  cf_assert(original.size() == 1);
  cf_assert(faces.size() == 1);

  ///Store T
  _transVector[0] = original[0];
  _coords[0] = coord;
  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();

  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

  ///Create the Geometric Entity to get the face
  ///corresponding to the idx/TRS given in .first and .second
  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  GeometricEntity& currFace = *geoBuilder.buildGE();

  /// Build the neighbor cell
  const CFuint faceID = currFace.getID();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;

  geoDataCell.trs = cellTrs;
  geoDataCell.idx = cellTrsID;
  GeometricEntity& neighborCell = *geoBuilderCell.buildGE();

  /// Compute the shape function gradient at the projected point
  const std::vector<RealMatrix> cellGradients = neighborCell.computeSolutionShapeFunctionGradients(_coords);
  const CFuint nbNodesInCell = neighborCell.getNodes()->size();

  /// Compute the gradient of T
  _gradientsT = 0.;
  for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode){
    for (CFuint iDim =0; iDim < _gradientsT.size() ; ++iDim){
      (_gradientsT)[iDim] += (cellGradients[0])(iNode,iDim) * (*(neighborCell.getState(iNode)))[0] ;
    }
  }

  /// Get the face normal
  ///@todo change this to use not the
  /// Average but the local FaceNormal (for 2nd order elements)
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  _normal = (neighborCell.computeAvgFaceNormals())[iFaceLocal];
  _normal.normalize();
/*const std::string namespaceName = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
CFout << "Namespace: " << namespaceName <<" - Pre Normal: " << _normal << "\n";*/
  /// Project the gradient on the face normal
  _transVector[1] = 0.;
  for (CFuint iDim = 0; iDim < _gradientsT.size() ; ++iDim){
    _transVector[1] += _gradientsT[iDim]*_normal[iDim];
  }
  _transVector[1] *= _model->getConductivity();


// std::cout << "-----------------------------------"<<std::endl;
// std::cout << "coord : " << coord <<std::endl;
// std::cout << "node0 coord : " << *(neighborCell.getNode(0)) <<std::endl;
// std::cout << "node1 coord : " << *(neighborCell.getNode(1)) <<std::endl;
// std::cout << "node2 coord : " << *(neighborCell.getNode(2)) <<std::endl;
// std::cout << "gradientsT[XX] : " << gradientsT[XX] <<std::endl;
// std::cout << "gradientsT[YY] : " << gradientsT[YY] <<std::endl;
// std::cout << "Normal : " << _normal <<std::endl;
// std::cout << "dTdn : " << _transVector[1] <<std::endl;

  ///release GeometricEntity
  geoBuilder.releaseGE();
  geoBuilderCell.releaseGE();

  return (&_transVector);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
